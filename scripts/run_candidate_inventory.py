#!/usr/bin/env python3
"""Freeze and inventory the full locked MIPLIB corpus with the C17 harness.

This is coverage exploration, never fair timing evidence. An explicit reviewed
commit and fresh output directory are required. Every locked case remains in
the derived inventory, including missing, crashed and malformed attempts.
"""
import argparse
import collections
import csv
import hashlib
import html
import json
import math
import os
import pathlib
import platform
import shutil
import statistics
import subprocess
import sys
import tarfile
import time

# A retained driver may be inspected/run beside its frozen evaluator helpers.
_helpers = pathlib.Path(__file__).resolve().parent / 'evaluator'
if (_helpers / 'benchmark_gate.py').is_file():
    sys.path.insert(0, str(_helpers))
from benchmark_gate import CORPUS_SIZE, validate_manifest, strict_json


def digest(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()


def write(path, value):
    with pathlib.Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def source(root):
    return {'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
            'dirty': bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=root, text=True))}


def unpack_source(archive, destination):
    """Extract the controlled git archive; reject links or path escapes."""
    with tarfile.open(archive) as stream:
        for member in stream.getmembers():
            parts = pathlib.PurePosixPath(member.name).parts
            if not parts or parts[0] != 'source' or '..' in parts or member.name.startswith('/'):
                raise ValueError('unexpected source archive path: ' + member.name)
            target = destination.joinpath(*parts)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open('xb') as output:
                    shutil.copyfileobj(stream.extractfile(member), output)
                target.chmod(member.mode & 0o777)
            else:
                raise ValueError('source archive links/special files are unsupported')
    return destination / 'source'


def tree_hashes(directory):
    return {str(path.relative_to(directory)): digest(path)
            for path in sorted(directory.rglob('*')) if path.is_file()}


def command(argv, directory, label, timeout):
    argv = list(map(str, argv))
    write(directory / (label + '.argv.json'), argv)
    started = time.monotonic()
    result = {'argv': argv, 'exit_code': None, 'timeout': False}
    with (directory / (label + '.stdout.txt')).open('w') as stdout, \
            (directory / (label + '.stderr.txt')).open('w') as stderr:
        try:
            result['exit_code'] = subprocess.run(argv, stdout=stdout, stderr=stderr,
                                                  timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            result['timeout'] = True
        except OSError as error:
            result['error'] = str(error)
    result['wall_s'] = time.monotonic() - started
    write(directory / (label + '.process.json'), result)
    return result


def summarize(manifest, records, logs=None, log_errors=None):
    """Pure accounting: missing and duplicate records cannot disappear."""
    logs = logs or {}
    log_errors = log_errors or {}
    grouped = collections.defaultdict(list)
    for record in records:
        grouped[record.get('instance')].append(record)
    locked = {model['id']: model for model in manifest['instances']}
    cases, missing, duplicates = [], [], []
    for name, model in locked.items():
        attempts = grouped.get(name, [])
        run = attempts[0] if len(attempts) == 1 else {}
        solver = logs.get(name) or {}
        checker = run.get('checker') or {}
        if not attempts:
            missing.append(name)
        elif len(attempts) != 1:
            duplicates.append(name)
        status = ('missing' if not attempts else 'duplicate' if len(attempts) != 1
                  else run.get('claimed_status', 'unknown'))
        cases.append({'instance': name, 'attempt_count': len(attempts),
                      'model_sha256': model['sha256'], 'observed_model_sha256': run.get('model_sha256'),
                      'reference_status': model['reference_status'],
                      'reference_objective': model.get('reference_objective'),
                      'status': status, 'message': solver.get('message'),
                      'exit_code': run.get('exit_code'), 'signal': run.get('signal'),
                      'reported_solver_status': run.get('reported_solver_status'),
                      'reported_solution_status': run.get('reported_solution_status'),
                      'status_conflict': run.get('status_conflict'),
                      'watchdog_timeout': run.get('watchdog_timeout'),
                      'solver_log_error': log_errors.get(name),
                      'rows': solver.get('rows'), 'columns': solver.get('columns'),
                      'nonzeros': solver.get('nonzeros'), 'nodes': solver.get('nodes'),
                      'optimizer_statistics': {key: solver[key] for key in
                          ['factor_core_peak', 'factor_nnz_peak', 'standard_nnz_peak',
                           'singleton_eliminations_peak', 'basis_updates', 'basis_update_rejections',
                           'presolve_passes', 'presolve_tightenings', 'presolve_fixed', 'presolve_rows']
                          if key in solver},
                      'checked_feasible': checker.get('status') == 'valid' and checker.get('checked_feasible') is True,
                      'checked_objective': checker.get('objective'), 'checker_status': checker.get('status'),
                      'run_record': f'runs/{name}/lab/seed-{run.get("seed", 1)}/run.json'})
    shape = {}
    for key in ['rows', 'columns', 'nonzeros']:
        values = [case[key] for case in cases if type(case[key]) is int]
        shape[key] = {'available': len(values), 'minimum': min(values) if values else None,
                      'median': statistics.median(values) if values else None,
                      'maximum': max(values) if values else None}
    return {'cases': cases, 'missing_instances': missing, 'duplicate_instances': duplicates,
            'unexpected_instances': sorted(str(name) for name in set(grouped)-set(locked)),
            'statuses': dict(collections.Counter(case['status'] for case in cases)),
            'stop_reasons': dict(collections.Counter(str(case['message']) for case in cases)),
            'shape': shape, 'checked_incumbents': sum(case['checked_feasible'] for case in cases),
            'malformed_or_missing_solver_logs': [{'instance': name, 'error': error}
                                                 for name, error in sorted(log_errors.items())],
            'invalid_candidates': sum((r.get('checker') or {}).get('status') == 'invalid' for r in records),
            'checker_errors': sum((r.get('checker') or {}).get('status') in ('error', 'unsupported') or
                                  r.get('checker_exit_code') not in (None, 0) and
                                  (r.get('checker') or {}).get('status') != 'invalid' for r in records),
            'crashes': [case['instance'] for case in cases if case['signal'] not in (None, 0) and
                        not case['watchdog_timeout']],
            'watchdog_limits': [case['instance'] for case in cases if case['watchdog_timeout']],
            'unknown_statuses': [case['instance'] for case in cases if case['status'] in ('unknown', None)],
            'model_hash_mismatches': [case['instance'] for case in cases if case['attempt_count'] == 1 and
                                      case['model_sha256'] != case['observed_model_sha256']]}


def html_report(path, summary):
    esc = lambda value: html.escape(str(value), quote=True)
    text = '''<!doctype html><meta charset="utf-8"><title>Full MIPLIB candidate inventory</title>
<style>body{font:16px system-ui;max-width:1450px;margin:36px auto;padding:0 24px;color:#123}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;padding:7px;border-bottom:1px solid #ccd}th{background:#123;color:white}a{color:#086d80}</style>
<h1>Full locked-corpus candidate inventory</h1>
<p>Coverage exploration with one solver process at a time. Concurrent engineering is permitted, so this is not a fair timing baseline, holdout or solver ranking. Solver statuses, checked witnesses, reference agreement and independent proofs remain distinct.</p>'''
    text += '<p>Frozen source <code>' + esc(summary['source_revision']) + '</code>; ' + \
            esc(summary['declared_limit_s']) + '-second cap, seed ' + esc(summary['seed']) + '.</p>'
    text += '<p>Retained membership: ' + esc(summary['membership_count']) + '; statuses: ' + esc(summary['statuses']) + \
            '. Checked witnesses: ' + esc(summary['checked_incumbents']) + '; correct completions: ' + \
            esc(summary.get('completed_correctly')) + '. Completion gate: ' + esc(summary.get('completion_gate_passed')) + '.</p>'
    text += '<p>Missing, duplicate, malformed, crashed, unsupported and limited cases remain unresolved. Unknown references are not upgraded to optima. No HiGHS speed wins are claimed; a future serial, quiescent paired campaign with matching budgets/hardware is required.</p>'
    text += '<p><a href="inventory.json">Inventory JSON</a> · <a href="instances.jsonl">Every retained case</a> · <a href="gate.json">Machine gate</a> · <a href="completion-manifest.json">Completion contracts</a> · <a href="inventory-provenance.json">Provenance</a> · <a href="evaluator-source.tar.gz">Frozen source</a> · <a href="runs/summary.html">Harness report</a></p>'
    text += '<table><tr><th>Case</th><th>Rows</th><th>Columns</th><th>Nonzeros</th><th>Status</th><th>Checked witness</th><th>Message</th><th>Evidence</th></tr>'
    for case in summary['cases']:
        text += '<tr>' + ''.join('<td>' + esc(case[key]) + '</td>' for key in
                               ['instance', 'rows', 'columns', 'nonzeros', 'status', 'checked_feasible', 'message'])
        text += '<td><a href="' + esc(case['run_record']) + '">record</a></td></tr>'
    with path.open('x') as stream:
        stream.write(text + '</table>')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--revision', required=True, help='reviewed committed source revision to freeze')
    parser.add_argument('--out', type=pathlib.Path, required=True)
    parser.add_argument('--limit', type=float, default=2)
    parser.add_argument('--grace', type=float, default=15)
    parser.add_argument('--seed', type=int, default=1)
    args = parser.parse_args()
    root = next(p for p in pathlib.Path(__file__).resolve().parents
                if (p / '.git').exists() and (p / 'CMakeLists.txt').exists())
    os.chdir(root)
    if not all(math.isfinite(x) and x > 0 for x in [args.limit, args.grace]) or args.seed <= 0:
        parser.error('positive finite time/grace and a positive seed are required')
    beginning = source(root)
    requested = subprocess.check_output(['git', 'rev-parse', args.revision + '^{commit}'], cwd=root, text=True).strip()
    if beginning['dirty'] or beginning['revision'] != requested:
        parser.error('freeze requires the requested committed revision and a clean working tree')
    locked_path = root / 'config/miplib.lock.json'
    manifest = strict_json(locked_path.read_text())
    _, errors = validate_manifest(manifest)
    if errors:
        parser.error('; '.join(errors))
    suite_path = root / 'data/miplib2017/suite.tsv'
    with suite_path.open() as stream:
        suite = list(csv.DictReader(stream, delimiter='\t'))
    locked = {r['id']: r for r in manifest['instances']}
    if len(suite) != CORPUS_SIZE or {r['id'] for r in suite} != set(locked):
        parser.error('suite membership must match the complete 240-model lock')
    for row in suite:
        expected = locked[row['id']]
        if row['sha256'] != expected['sha256'] or digest(root / row['path']) != expected['sha256']:
            parser.error('dataset hash mismatch for ' + row['id'])
        if row['reference_status'] != expected['reference_status'] or row['reference_objective'] != expected['reference_objective']:
            parser.error('suite reference differs from the dataset lock for ' + row['id'])
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    subprocess.run(['git', 'archive', '--format=tar.gz', '--prefix=source/', requested,
                    '-o', str(out / 'evaluator-source.tar.gz')], cwd=root, check=True)
    frozen_source = unpack_source(out / 'evaluator-source.tar.gz', out)
    initial_source_hashes = tree_hashes(frozen_source)
    write(out / 'archived-source-hashes.json', initial_source_hashes)
    evaluator = out / 'evaluator'
    evaluator.mkdir()
    config = out / 'config'
    shutil.copytree(frozen_source / 'config', config)
    if digest(config / 'miplib.lock.json') != digest(locked_path):
        raise RuntimeError('dataset lock changed between preflight and archive')
    shutil.copy2(suite_path, out / 'suite.tsv')
    shutil.copy2(__file__, out / 'driver.py')
    build = out / 'build'
    configure = command(['cmake', '-S', frozen_source, '-B', build, '-DCMAKE_BUILD_TYPE=Release'],
                        out, 'configure', 120)
    compiled = command(['cmake', '--build', build, '--target', 'milpsolve', 'milpcheck', 'milpbench',
                        '--parallel', '4'], out, 'build', 300) if configure['exit_code'] == 0 else None
    if configure['exit_code'] != 0 or not compiled or compiled['exit_code'] != 0:
        write(out / 'pipeline-failure.json', {'phase': 'building archived evaluator',
              'source_revision': requested, 'configure': configure, 'build': compiled,
              'solver_processes_launched': 0, 'unattempted_instances': sorted(locked)})
        return 1
    cache = (build / 'CMakeCache.txt').read_text()
    compiler = next(line.split('=', 1)[1] for line in cache.splitlines()
                    if line.startswith('CMAKE_C_COMPILER:'))
    compiler_identity = command([compiler, '--version'], out, 'compiler-version', 30)
    identities = {}
    for name in ['milpsolve', 'milpcheck', 'milpbench']:
        original = build / name
        target = evaluator / name
        shutil.copy2(original, target)
        identities[name] = {'origin': str(original), 'path': str(target), 'sha256': digest(target)}
        if digest(original) != digest(target):
            raise RuntimeError('binary changed while freezing: ' + name)
    for name in ['benchmark_gate.py', 'report.py']:
        shutil.copy2(frozen_source / 'scripts' / name, evaluator / name)
    (out / 'solvers.tsv').write_text('solver\tbinary\tsha256\tversion\nlab\t' + str(evaluator / 'milpsolve') +
        '\t' + identities['milpsolve']['sha256'] + '\t' + requested[:12] + ' coverage inventory\n')
    manifest['validation_contract_version'] = 'miplib-reference-completion-v1'
    manifest['derived_from_dataset_lock_sha256'] = digest(locked_path)
    for model in manifest['instances']:
        if model['reference_status'] == '=inf=' and 'validation_contract' not in model:
            model['validation_contract'] = 'reference_status_agreement'
            model['validation_basis'] = {'reference': manifest.get('reference'),
                'scope': 'Declared benchmark expected-status agreement; not an independent proof.'}
    write(out / 'completion-manifest.json', manifest)
    provenance = {'schema': 1, 'kind': 'Full locked-corpus development coverage inventory; not fair timing or ranking',
                  'coverage_exploration': True, 'fair_timing_campaign': False, 'source_begin': beginning,
                  'source_revision': requested, 'source_archive_sha256': digest(out / 'evaluator-source.tar.gz'),
                  'driver_sha256': digest(out / 'driver.py'), 'identities': identities,
                  'evaluators_built_from_archive': True, 'configure': configure, 'build': compiled,
                  'cmake_cache_sha256': digest(build / 'CMakeCache.txt'),
                  'compiler_identity': compiler_identity,
                  'compiler_version_output_sha256': digest(out / 'compiler-version.stdout.txt'),
                  'archived_source_hashes_sha256': digest(out / 'archived-source-hashes.json'),
                  'membership_count': CORPUS_SIZE, 'time_limit_s': args.limit, 'grace_s': args.grace,
                  'seed': args.seed, 'requested_threads': 1, 'serial_solver_processes': True,
                  'suite_sha256': digest(out / 'suite.tsv'), 'solver_config_sha256': digest(out / 'solvers.tsv'),
                  'configuration_hashes': {p.name: digest(p) for p in config.iterdir() if p.is_file()},
                  'host': {'system': platform.system(), 'machine': platform.machine(), 'release': platform.release()}}
    write(out / 'inventory-provenance.json', provenance)
    print('Frozen', requested, 'all', CORPUS_SIZE, 'models; starting serial coverage inventory.', flush=True)
    run = command([evaluator / 'milpbench', '--mip', '--suite', out / 'suite.tsv',
                   '--solvers', out / 'solvers.tsv', '--checker', evaluator / 'milpcheck',
                   '--limit', str(args.limit), '--grace', str(args.grace), '--seed', str(args.seed),
                   '--out', out / 'runs'], out, 'runner', CORPUS_SIZE*(args.limit+args.grace+2)+30)
    raw = out / 'runs/runs.jsonl'
    records, malformed = [], []
    if raw.exists():
        for line_number, line in enumerate(raw.read_text().splitlines(), 1):
            if not line.strip():
                continue
            try:
                value = strict_json(line)
                if not isinstance(value, dict):
                    raise ValueError('run record is not an object')
                if not isinstance(value.get('instance'), str) or not isinstance(value.get('checker'), dict):
                    raise ValueError('run record instance/checker schema is malformed')
                records.append(value)
            except ValueError as error:
                malformed.append({'line': line_number, 'error': str(error)})
    logs, log_errors = {}, {}
    for row in records:
        path = out / 'runs' / str(row.get('instance')) / 'lab' / f'seed-{args.seed}' / 'solver.log'
        try:
            value = strict_json(path.read_text())
            if not isinstance(value, dict) or not isinstance(value.get('status'), str):
                raise ValueError('solver log is not a solver status JSON object')
            logs[row.get('instance')] = value
        except (OSError, ValueError) as error:
            log_errors[row.get('instance')] = str(error)
    summary = summarize(manifest, records, logs, log_errors)
    if raw.exists():
        gate_process = command([sys.executable, evaluator / 'benchmark_gate.py', '--manifest',
            out / 'completion-manifest.json', '--runs', raw, '--out', out / 'gate.json'], out, 'gate', 120)
    else:
        gate_process = {'exit_code': None, 'error': 'run records missing'}
        write(out / 'gate.json', {'passed': False, 'completion_gate_passed': False,
                                  'input_error': 'run records missing'})
    try:
        gate = strict_json((out / 'gate.json').read_text())
        if not isinstance(gate, dict):
            raise ValueError('gate did not emit a JSON object')
    except (OSError, ValueError) as error:
        gate = {'passed': False, 'completion_gate_passed': False,
                'input_error': 'gate artifact unavailable: ' + str(error)}
        write(out / ('gate-failure.json' if (out / 'gate.json').exists() else 'gate.json'), gate)
    harness_report = None
    if records and not malformed:
        command([sys.executable, evaluator / 'report.py', out / 'runs'], out, 'harness-report', 120)
        try:
            harness_report = strict_json((out / 'runs/summary.json').read_text())
        except (OSError, ValueError):
            pass
    end_hashes = {name: digest(pathlib.Path(identity['path'])) for name, identity in identities.items()}
    summary.update(schema=1, source_revision=requested, coverage_exploration=True,
        fair_timing_campaign=False, membership_count=CORPUS_SIZE, attempt_records=len(records),
        declared_limit_s=args.limit, seed=args.seed, malformed_run_lines=malformed, runner=run,
        gate_error=gate.get('input_error'),
        gate_process=gate_process, completed_correctly=gate.get('completed_correctly'),
        completion_gate_passed=gate.get('completion_gate_passed', False), speed_gate_passed=False,
        strict_candidate_wins=0, reference_status_counts=dict(collections.Counter(r['reference_status'] for r in manifest['instances'])),
        unresolved_reference_instances=gate.get('unresolved_reference_instances', []),
        source_end=source(root), frozen_binary_hashes_end=end_hashes)
    summary['archived_source_unchanged'] = tree_hashes(frozen_source) == initial_source_hashes
    summary['reference_contradictions'] = ((harness_report or {}).get('summary', {}).get('lab', {})
                                           .get('reference_contradictions'))
    summary['inventory_complete'] = (len(records) == CORPUS_SIZE and not malformed and
        not summary['missing_instances'] and not summary['duplicate_instances'] and
        not summary['unexpected_instances'] and not summary['model_hash_mismatches'] and
        run['exit_code'] == 0 and summary['archived_source_unchanged'] and
        all(end_hashes[name] == identity['sha256'] for name, identity in identities.items()))
    summary['largest_by_nonzeros'] = sorted((r for r in summary['cases'] if type(r['nonzeros']) is int),
                                           key=lambda r: r['nonzeros'], reverse=True)[:20]
    summary['smallest_unsupported_by_original_rows_plus_columns'] = sorted(
        (r for r in summary['cases'] if r['status'] == 'unsupported' and
         type(r['rows']) is int and type(r['columns']) is int), key=lambda r: r['rows']+r['columns'])[:20]
    write(out / 'inventory.json', summary)
    with (out / 'instances.jsonl').open('x') as stream:
        for case in summary['cases']:
            stream.write(json.dumps(case) + '\n')
    html_report(out / 'report.html', summary)
    write(out / 'inventory-final-lock.json', {'source_revision': requested, 'frozen_binary_hashes': end_hashes,
        'artifacts': {name: digest(out / name) for name in ['inventory.json', 'instances.jsonl', 'gate.json',
                     'completion-manifest.json', 'inventory-provenance.json', 'suite.tsv', 'solvers.tsv', 'driver.py']}})
    print(json.dumps({key: summary[key] for key in ['inventory_complete', 'statuses', 'checked_incumbents',
        'invalid_candidates', 'checker_errors', 'completed_correctly', 'completion_gate_passed']}), flush=True)
    return 0 if summary['inventory_complete'] else 1


if __name__ == '__main__':
    sys.exit(main())
