#!/usr/bin/env python3
"""Controlled root-LP development comparison for scaling and refinement.

The three models were selected because previous development runs exposed
numerical failures. They are not a holdout. Every ablation and failure is
retained; only the default new engine and the HiGHS reference are acceptance
gates. Cross-solver objective agreement is not an optimality certificate.
"""
import argparse
import collections
import csv
import datetime
import hashlib
import json
import math
import os
import pathlib
import platform
import re
import shutil
import signal
import subprocess
import time
from fractions import Fraction as F

ROOT = pathlib.Path(__file__).resolve().parents[1]
MODELS = ('mas74', 'mas76', 'gen-ip054')
PROFILES = (
    ('old', 'old', ()),
    ('unscaled_unrefined', 'new', ('--no-scaling', '--no-refinement')),
    ('scaling_only', 'new', ('--no-refinement',)),
    ('refinement_only', 'new', ('--no-scaling',)),
    ('both', 'new', ()),
)


def absolute(path):
    path = pathlib.Path(path)
    return (path if path.is_absolute() else ROOT / path).resolve()


def digest(path):
    value = hashlib.sha256()
    with pathlib.Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(part)
    return value.hexdigest()


def source_identity():
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT,
                                        text=True).strip())
    paths = sorted([ROOT / 'CMakeLists.txt', *ROOT.glob('src/*.c'), *ROOT.glob('include/*.h'),
                    ROOT / 'config/miplib.lock.json', ROOT / 'config/solvers.tsv',
                    pathlib.Path(__file__).resolve()])
    return {'git_revision': revision, 'working_tree_dirty': dirty,
            'files': {str(path.relative_to(ROOT)): digest(path) for path in paths}}


def reject_nonfinite(value):
    raise ValueError('nonfinite JSON number: ' + value)


def read_json(text):
    try:
        value = json.loads(text, parse_constant=reject_nonfinite)
    except (ValueError, TypeError):
        return None
    def finite(item):
        if isinstance(item, float):
            return math.isfinite(item)
        if isinstance(item, dict):
            return all(finite(child) for child in item.values())
        if isinstance(item, list):
            return all(finite(child) for child in item)
        return True
    return value if isinstance(value, dict) and finite(value) else None


def run_process(command, folder, name, timeout):
    """Run serially with a process-group watchdog; preserve stdout and stderr."""
    environment = os.environ.copy()
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        environment[key] = '1'
    started = time.monotonic()
    killed = False
    launch_error = None
    stdout = stderr = ''
    try:
        child = subprocess.Popen(command, cwd=ROOT, env=environment, text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 start_new_session=True)
    except OSError as error:
        launch_error = str(error)
        code = None
    else:
        try:
            stdout, stderr = child.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            killed = True
            try:
                os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                stdout, stderr = child.communicate(timeout=1)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                stdout, stderr = child.communicate()
        code = child.returncode
    (folder / (name + '.stdout')).write_text(stdout)
    (folder / (name + '.stderr')).write_text(stderr)
    return {'argv': command, 'exit_code': code,
            'signal': -code if code is not None and code < 0 else 0,
            'watchdog_timeout': killed, 'launch_error': launch_error,
            'wall_s': time.monotonic() - started,
            'stdout': name + '.stdout', 'stderr': name + '.stderr'}, stdout, stderr


def highs_status(stdout):
    status = 'unknown'
    for line in stdout.splitlines():
        match = re.match(r'^\s*Model status\s*:\s*(.*?)\s*$', line)
        if match:
            status = {'Optimal': 'optimal', 'Infeasible': 'infeasible',
                      'Unbounded': 'unbounded', 'Time limit reached': 'time_limit',
                      'Infeasible or unbounded': 'infeasible_or_unbounded'}.get(match[1], 'unknown')
    return status


def check_candidate(checker, model, folder, fmt):
    solution, dual = folder / 'solution.sol', folder / 'dual.sol'
    if not solution.exists() and not dual.exists():
        return {'status': 'unavailable', 'checked_feasible': False,
                'bound_status': 'absent', 'certifies_optimality': False}, None
    command = [str(checker), '--model', str(model), '--relax']
    if solution.exists():
        command += ['--solution', str(solution), '--format', fmt]
    if dual.exists():
        command += ['--dual', str(dual)]
    process, stdout, _ = run_process(command, folder, 'checker', 120)
    result = read_json(stdout)
    if result is None:
        result = {'status': 'invalid_output', 'checked_feasible': False,
                  'bound_status': 'absent', 'certifies_optimality': False}
    (folder / 'checker.json').write_text(json.dumps(result, indent=2) + '\n')
    return result, process


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', required=True, help='fresh artifact directory')
    ap.add_argument('--solver', default='build/milpsolve')
    ap.add_argument('--old-solver', default='results/mip-native-v2/evaluator/milpsolve')
    ap.add_argument('--checker', default='build/milpcheck')
    ap.add_argument('--highs', default='tools/highs/bin/highs')
    ap.add_argument('--limit', type=float, default=5)
    ap.add_argument('--iteration-limit', type=int, default=5000)
    args = ap.parse_args()
    if not math.isfinite(args.limit) or args.limit <= 0 or args.iteration_limit <= 0:
        ap.error('positive finite time and iteration limits are required')
    out = absolute(args.out)
    if out.exists():
        ap.error('--out must not exist; previous experiments are immutable')
    current, old, checker, highs = map(absolute,
                                      (args.solver, args.old_solver, args.checker, args.highs))
    for path in (current, old, checker, highs):
        if not path.is_file():
            ap.error('required input is missing: ' + str(path))
    if digest(current) == digest(old):
        ap.error('old and new solvers must be distinct executable snapshots')
    with (ROOT / 'config/solvers.tsv').open() as stream:
        high_config = next(row for row in csv.DictReader(stream, delimiter='\t')
                           if row['solver'] == 'highs')
    high_sha = digest(highs)
    if high_sha != high_config['sha256']:
        ap.error('HiGHS executable does not match the pinned solver digest')
    manifest_path = ROOT / 'config/miplib.lock.json'
    manifest = json.loads(manifest_path.read_text())
    locked = {row['id']: row for row in manifest['instances']}
    models = {}
    for name in MODELS:
        path = ROOT / 'data/miplib2017/instances' / (name + '.mps.gz')
        if not path.is_file() or digest(path) != locked[name]['sha256']:
            ap.error('missing or changed locked model: ' + name)
        models[name] = {'path': str(path), 'sha256': locked[name]['sha256']}
    beginning = source_identity()
    out.mkdir(parents=True)
    evaluator = out / 'evaluator'; evaluator.mkdir()
    origins = {'new': str(current), 'old': str(old), 'checker': str(checker)}
    for path, name in ((current, 'milpsolve-new'), (old, 'milpsolve-old'), (checker, 'milpcheck')):
        shutil.copy2(path, evaluator / name)
    binaries = {'new': evaluator / 'milpsolve-new', 'old': evaluator / 'milpsolve-old'}
    checker = evaluator / 'milpcheck'
    identities = {'origins': origins,
                  'frozen': {name: {'path': str(path), 'sha256': digest(path)}
                             for name, path in {**binaries, 'checker': checker}.items()},
                  'highs': {'path': str(highs), 'sha256': high_sha, 'version': high_config['version'],
                            'copied': False, 'reason': 'preserve the release installation and library lookup paths'}}
    old_manifest = old.parent.parent / 'summary.json'
    identities['old_recorded_source'] = None
    if old_manifest.is_file():
        prior = read_json(old_manifest.read_text())
        if prior and prior.get('binaries', {}).get('milpsolve') == identities['frozen']['old']['sha256']:
            identities['old_recorded_source'] = {'manifest': str(old_manifest),
                'manifest_sha256': digest(old_manifest), 'source': prior.get('source')}
    records = []
    for name in MODELS:
        model = pathlib.Path(models[name]['path'])
        if digest(model) != models[name]['sha256']:
            raise RuntimeError('model changed during experiment preparation: ' + name)
        for profile, executable, flags in PROFILES:
            folder = out / name / profile; folder.mkdir(parents=True)
            command = [str(binaries[executable]), str(model), '--relax', '--time-limit', str(args.limit),
                       '--iteration-limit', str(args.iteration_limit), '--solution',
                       str(folder / 'solution.sol'), '--dual', str(folder / 'dual.sol'), *flags]
            (folder / 'settings.json').write_text(json.dumps({'profile': profile, 'relaxation': True,
                'time_limit_s': args.limit, 'iteration_limit': args.iteration_limit,
                'requested_threads': 1, 'flags': list(flags)}, indent=2) + '\n')
            process, stdout, _ = run_process(command, folder, 'solver', args.limit + 15)
            result = read_json(stdout)
            status = result.get('status', 'invalid_output') if result else 'invalid_output'
            check, check_process = check_candidate(checker, model, folder, 'simple')
            record = {'instance': name, 'profile': profile, 'status': status,
                      'process': process, 'solver': result, 'checker': check,
                      'checker_process': check_process, 'model': models[name],
                      'binary_sha256': identities['frozen'][executable]['sha256'],
                      'record': str((folder / 'run.json').relative_to(out))}
            (folder / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
            records.append(record)
            print(f'{name} {profile}: {status}; checker={check.get("status")}', flush=True)
        folder = out / name / 'highs'; folder.mkdir(parents=True)
        options = folder / 'options.txt'
        options.write_text('\n'.join([
            'solve_relaxation = true', 'threads = 1', 'random_seed = 1',
            f'time_limit = {args.limit}', 'write_solution_to_file = true',
            'write_solution_style = 0', f'solution_file = {folder / "solution.sol"}',
            f'log_file = {folder / "highs.log"}', '']) )
        command = [str(highs), str(model), '--options_file', str(options)]
        process, stdout, _ = run_process(command, folder, 'solver', args.limit + 15)
        status = highs_status(stdout)
        check, check_process = check_candidate(checker, model, folder, 'highs')
        record = {'instance': name, 'profile': 'highs', 'status': status, 'process': process,
                  'solver': None, 'checker': check, 'checker_process': check_process,
                  'model': models[name], 'binary_sha256': high_sha,
                  'record': str((folder / 'run.json').relative_to(out))}
        (folder / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
        records.append(record)
        print(f'{name} highs: {status}; checker={check.get("status")}', flush=True)
    failures, comparisons = [], []
    by_pair = {(r['instance'], r['profile']): r for r in records}
    if len(records) != 18 or len(by_pair) != 18:
        failures.append('missing, duplicate or unexpected profile/model records')
    for name in MODELS:
        reference, candidate = by_pair[name, 'highs'], by_pair[name, 'both']
        for record in (reference, candidate):
            label = name + '/' + record['profile']
            process, check = record['process'], record['checker']
            if process['exit_code'] != 0 or process['signal'] or process['watchdog_timeout']:
                failures.append(label + ': process did not exit cleanly')
            if record['status'] != 'optimal':
                failures.append(label + ': solver status ' + record['status'])
            if not (check.get('status') == 'valid' and check.get('checked_feasible')):
                failures.append(label + ': original-model primal checking did not pass')
            if record['checker_process'] is None or record['checker_process']['exit_code'] != 0:
                failures.append(label + ': checker process did not exit cleanly')
            if record['solver'] and record['solver'].get('certifies_optimality', False):
                failures.append(label + ': optimizer made an unsupported independent-proof claim')
        a, b = candidate['checker'].get('objective'), reference['checker'].get('objective')
        difference = abs(F(a) - F(b)) if a is not None and b is not None else None
        tolerance = F(1, 10**6) + F(1, 10**9) * abs(F(b)) if b is not None else None
        agrees = difference is not None and difference <= tolerance
        if not agrees:
            failures.append(name + ': default new objective differs from checked HiGHS root objective')
        comparisons.append({'instance': name, 'new_checked_objective': a,
                            'highs_checked_objective': b, 'absolute_difference': str(difference)
                            if difference is not None else None,
                            'comparison_tolerance': str(tolerance) if tolerance is not None else None,
                            'objectives_agree': agrees,
                            'new_bound_status': candidate['checker'].get('bound_status'),
                            'new_certified_bound': candidate['checker'].get('certified_bound'),
                            'new_certifies_optimality': candidate['checker'].get('certifies_optimality', False)})
    for name, model in models.items():
        if digest(model['path']) != model['sha256']:
            failures.append('model changed during experiment: ' + name)
    if digest(highs) != high_sha:
        failures.append('pinned HiGHS executable changed during experiment')
    ending = source_identity()
    if beginning != ending:
        failures.append('source revision or tracked source files changed during experiment')
    if beginning['working_tree_dirty'] or ending['working_tree_dirty']:
        failures.append('experiment source was not a clean Git revision')
    summary = {'schema': 1, 'kind': 'declared development root-LP numerical comparison; not a speed ranking',
               'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'selection': 'mas74/mas76 basis-residual failures and gen-ip054 Phase-I development failure',
               'holdout': False, 'models': models, 'dataset_manifest_sha256': digest(manifest_path),
               'profiles': [{'name': name, 'executable': executable, 'flags': list(flags)}
                            for name, executable, flags in PROFILES] + [{'name': 'highs', 'relaxation': True}],
               'internal_time_limit_s': args.limit, 'iteration_limit': args.iteration_limit,
               'watchdog_policy': 'time limit + 15 seconds; SIGTERM then SIGKILL after one second',
               'requested_threads': 1, 'seed': 1,
               'objective_comparison': 'absolute 1e-6 + relative 1e-9 against independently checked HiGHS root solution',
               'proof_scope': 'Cross-solver agreement is not proof. Independent rational LP-bound status is separate; primal checks use recorded tolerances.',
               'passed': not failures, 'failures': failures, 'expected_runs': 18, 'actual_runs': len(records),
               'identities': identities, 'source_begin': beginning, 'source_end': ending,
               'platform': {'system': platform.system(), 'machine': platform.machine(), 'release': platform.release()},
               'counts': {'statuses': dict(collections.Counter(r['status'] for r in records)),
                          'checked_primal': sum(bool(r['checker'].get('checked_feasible')) for r in records),
                          'certified_lp_bounds': sum(r['checker'].get('bound_status') == 'valid' for r in records)},
               'comparisons': comparisons,
               'records': [{'instance': r['instance'], 'profile': r['profile'], 'status': r['status'],
                            'exit_code': r['process']['exit_code'], 'watchdog_timeout': r['process']['watchdog_timeout'],
                            'checked_feasible': bool(r['checker'].get('checked_feasible')),
                            'checked_objective': r['checker'].get('objective'),
                            'bound_status': r['checker'].get('bound_status'),
                            'certified_bound': r['checker'].get('certified_bound'),
                            'solver_diagnostics': r['solver'], 'record': r['record']}
                           for r in records]}
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({key: summary[key] for key in ('passed', 'actual_runs', 'counts', 'failures')}, indent=2))
    return 0 if summary['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
