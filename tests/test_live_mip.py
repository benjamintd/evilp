#!/usr/bin/env python3
"""Optional retained correctness crosscheck: native MILP baselines and lab B&B.

All fixture optima are independently computed by exact enumeration. This is
not a performance ranking; raw artifacts, exits and all failures are retained.
"""
import argparse
import csv
import datetime
import hashlib
import json
import pathlib
import shutil
import subprocess
from fractions import Fraction as F

from test_mip_solver import integer_mps, mip_oracle

ROOT = pathlib.Path(__file__).resolve().parents[1]


def sha(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(part)
    return digest.hexdigest()


def fixtures():
    return [
        ('binary_gap', {'cost': [1, 1], 'integer': [True, True],
                        'bounds': [(0, 1), (0, 1)], 'rows': [([2, 2], None, 3)], 'maximize': True}),
        ('general_integer', {'cost': [-2, 3], 'integer': [True, True],
                             'bounds': [(-2, 2), (-1, 2)], 'rows': [([1, 1], F(1, 2), None)]}),
        ('mixed_fraction', {'cost': [-2, 1], 'integer': [True, False],
                            'bounds': [(0, 3), (0, 2)], 'rows': [([1, 3], 4, None)], 'offset': F(7, 4)}),
        ('max_offset_fixed', {'cost': [3, 2, F(1, 10)], 'integer': [True, True, False],
                              'bounds': [(0, 1), (0, 1), (F(3, 10), F(3, 10))],
                              'rows': [([2, 1, 0], None, 2)], 'maximize': True, 'offset': 5}),
        ('fractional_bounds', {'cost': [-1, 2], 'integer': [True, True],
                               'bounds': [(F(-3, 2), F(5, 2)), (-1, 2)],
                               'rows': [([1, 1], -1, F(5, 2))]}),
        ('mixed_free', {'cost': [2, 1], 'integer': [True, False],
                        'bounds': [(-2, 2), (None, None)],
                        'rows': [([1, 3], 1, 1), ([0, 1], -2, 2)], 'maximize': True}),
        ('fixed_redundant', {'cost': [1, 2, -1], 'integer': [True, True, False],
                             'bounds': [(0, 3), (1, 1), (0, 3)],
                             'rows': [([1, 0, 1], 2, 2), ([2, 0, 2], 4, 4)]}),
        ('integer_infeasible', {'cost': [1], 'integer': [True],
                                'bounds': [(0, 1)], 'rows': [([1], F(1, 2), F(1, 2))]}),
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', required=True)
    ap.add_argument('--runner', default='build/milpbench')
    ap.add_argument('--checker', default='build/milpcheck')
    ap.add_argument('--solver', default='build/milpsolve')
    ap.add_argument('--native-solvers', default='config/solvers.tsv')
    args = ap.parse_args()
    absolute = lambda p: (ROOT / p).resolve()
    out = absolute(args.out)
    if out.exists():
        ap.error('--out must be a fresh directory')
    runner, checker, solver = map(absolute, (args.runner, args.checker, args.solver))
    config = absolute(args.native_solvers)
    for path in (runner, checker, solver, config):
        if not path.is_file():
            ap.error('missing required input: ' + str(path))
    with config.open() as stream:
        native = list(csv.DictReader(stream, delimiter='\t'))
    if len(native) != 3 or {row['solver'] for row in native} != {'cbc', 'scip', 'highs'}:
        ap.error('native config must contain exactly HiGHS, SCIP and Cbc')
    models = [(name, model, mip_oracle(model)) for name, model in fixtures()]
    out.mkdir(parents=True)
    evaluator = out / 'evaluator'; evaluator.mkdir()
    for path, name in ((runner, 'milpbench'), (checker, 'milpcheck'), (solver, 'milpsolve')):
        shutil.copy2(path, evaluator / name)
    runner, checker, solver = [evaluator / x for x in ('milpbench', 'milpcheck', 'milpsolve')]
    model_dir = out / 'models'; model_dir.mkdir()
    expected = {}
    with (out / 'suite.tsv').open('x') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(['id', 'path', 'sha256', 'reference_status', 'reference_objective'])
        for name, model, value in models:
            path = model_dir / (name + '.mps'); path.write_text(integer_mps(model))
            expected[name] = {'objective': str(value) if value is not None else None,
                              'maximize': bool(model.get('maximize')), 'model_sha256': sha(path)}
            writer.writerow([name, path, sha(path), '=inf=' if value is None else '=opt=',
                             '-' if value is None else str(value)])
    shutil.copy2(config, out / 'native-solvers.tsv')
    command = [str(runner), '--suite', str(out / 'suite.tsv'), '--solvers',
               str(out / 'native-solvers.tsv'), '--checker', str(checker),
               '--out', str(out / 'native'), '--limit', '10', '--seed', '1']
    proc = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, timeout=300)
    (out / 'native-driver.log').write_text(proc.stdout + proc.stderr)
    failures, rows = [], []
    if proc.returncode:
        failures.append(f'native harness exit {proc.returncode}')
    native_runs = out / 'native/runs.jsonl'
    if native_runs.exists():
        rows = [json.loads(line) for line in native_runs.read_text().splitlines()]
    else:
        failures.append('missing native run records')
    for name, model, reference in models:
        for cold in (False, True):
            label = 'lab-cold' if cold else 'lab'
            folder = out / label / name; folder.mkdir(parents=True)
            solution = folder / 'solution.sol'
            command = [str(solver), str(model_dir / (name + '.mps')), '--mip',
                       '--solution', str(solution), '--time-limit', '10'] + (['--cold'] if cold else [])
            p = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, timeout=30)
            (folder / 'solver.json').write_text(p.stdout)
            (folder / 'solver.stderr').write_text(p.stderr)
            try:
                result = json.loads(p.stdout)
            except ValueError:
                failures.append(f'{name}/{label}: non-JSON solver output')
                continue
            check = {'checked_feasible': False, 'status': 'no_solution'}
            if solution.exists():
                q = subprocess.run([str(checker), '--model', str(model_dir / (name + '.mps')),
                                    '--solution', str(solution)], text=True, capture_output=True, timeout=15)
                (folder / 'checker.json').write_text(q.stdout)
                try:
                    check = json.loads(q.stdout)
                except ValueError:
                    failures.append(f'{name}/{label}: non-JSON checker output')
            record = {'instance': name, 'solver': label, 'claimed_status': result.get('status'),
                      'exit_code': p.returncode, 'checker': check, 'argv': command,
                      'model_sha256': expected[name]['model_sha256'], 'solver_sha256': sha(solver),
                      'checker_sha256': sha(checker), 'optimizer': result,
                      'record': str((folder / 'run.json').relative_to(out))}
            (folder / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
            rows.append(record)
    comparisons = []
    for row in rows:
        name, engine = row['instance'], row['solver']
        reference = expected[name]['objective']
        correct_status = row['claimed_status'] == ('infeasible' if reference is None else 'optimal')
        check = row['checker']
        feasible = bool(check.get('checked_feasible') and check.get('status') == 'valid')
        objective = F(check['objective']) if check.get('objective') is not None else None
        objective_matches = (not feasible if reference is None else objective is not None and
                             abs(objective - F(reference)) <= F(1, 10**7) + F(1, 10**9) * abs(F(reference)))
        if not correct_status or not objective_matches or row['exit_code'] != 0:
            failures.append(f'{name}/{engine}: status={row["claimed_status"]},exit={row["exit_code"]},objective={objective},reference={reference}')
        if reference is not None and not feasible:
            failures.append(f'{name}/{engine}: original-model candidate fails checker')
        if engine.startswith('lab') and (row['optimizer'].get('certifies_optimality') or
                                          row['optimizer'].get('bound_certified')):
            failures.append(f'{name}/{engine}: unsupported independent-proof claim')
        comparisons.append({'instance': name, 'solver': engine, 'claimed_status': row['claimed_status'],
                            'exit_code': row['exit_code'], 'passed_checker': feasible,
                            'objective': str(objective) if objective is not None else None,
                            'reference': reference, 'objective_matches': objective_matches,
                            'reference_outcome_matches': bool(correct_status and objective_matches)})
    pairs = [(r['instance'], r['solver']) for r in rows]
    wanted = {(name, engine) for name, _, _ in models for engine in ('cbc', 'scip', 'highs', 'lab', 'lab-cold')}
    if len(pairs) != len(set(pairs)) or set(pairs) != wanted:
        failures.append('missing, duplicate or unexpected solver/model pairs')
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True).strip())
    summary = {'schema': 1, 'kind': 'tiny native MILP correctness crosscheck; not a speed ranking',
               'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'passed': not failures, 'failures': failures, 'instances': len(models),
               'expected_runs': len(wanted), 'actual_runs': len(rows), 'reference_models': expected,
               'counts': {'matched_reference_outcomes': sum(r['reference_outcome_matches'] for r in comparisons),
                          'checked_feasible_candidates': sum(r['passed_checker'] for r in comparisons)},
               'source': {'git_revision': revision, 'working_tree_dirty': dirty,
                          'script_sha256': sha(__file__), 'fixture_helper_sha256': sha(ROOT / 'tests/test_mip_solver.py')},
               'binaries': {p.name: sha(p) for p in (runner, checker, solver)}, 'records': comparisons}
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: summary[k] for k in ('passed', 'instances', 'actual_runs', 'counts', 'failures')}, indent=2))
    return 0 if summary['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
