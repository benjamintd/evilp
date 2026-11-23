"""Run the six locked private fleet snapshots and replay the archived false claims."""
import argparse
import csv
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys

from run_evaluation import assess_evidence, checker_policy_errors, read_witness_ledgers, strict_json, witness_index
from verify_candidate import ROOT, TOOLS, sha, verified


def write(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def one(rows, **fields):
    matches = [row for row in rows if all(row.get(k) == v for k, v in fields.items())]
    if len(matches) != 1:
        raise ValueError(f'expected one record matching {fields}, found {len(matches)}')
    return matches[0]


def checked_point(row, model, path, fmt, checker):
    if sha(path) != row['evidence_sha256']['solution.sol']:
        raise ValueError(f'archived solution changed: {path}')
    cmd = [checker, '--model', model['path'], '--solution', path, '--format', fmt]
    observed = strict_json(subprocess.check_output(cmd))
    if (not observed.get('checked_feasible') or observed.get('relaxation') is not False
            or checker_policy_errors(observed) or observed['objective'] != row['checker']['objective']):
        raise ValueError(f'archived witness failed original-model replay: {path}')
    return {**row, 'checker': observed, 'checker_sha256': sha(checker), 'checker_exit_code': 0}


def replay_controls(config, models, checker, archived):
    pair_file = ROOT / 'results/fleet-repair-60s-v1/runs/runs.jsonl'
    pairs, _ = read_witness_ledgers([pair_file])
    panel = ROOT / config['legacy_panel']
    if sha(panel) != config['legacy_panel_sha256']:
        raise ValueError('historical panel changed')
    old = one([strict_json(line) for line in panel.read_text().splitlines()], id='fleet-2330')
    controls = []
    for suffix, better_solver in [('2300', 'highs'), ('2330', 'lab')]:
        name = 'fleet-20261117T' + suffix
        model = models[name]
        better = one(pairs, instance=name, solver=better_solver)
        def paired_point(row):
            path = pair_file.parent / name / row['solver'] / 'seed-1/solution.sol'
            return checked_point(row, model, path, 'highs' if row['solver'] == 'highs' else 'simple', checker)
        better = paired_point(better)
        expected = {'optimal_objective_worse_than_checked_witness'}
        if suffix == '2300':
            claim = paired_point(one(pairs, instance=name, solver='lab'))
        else:
            path = panel.parent / 'solutions/fleet-2330.sol'
            row = one(archived, instance=name, solution_path=str(path.relative_to(ROOT)))
            claim = checked_point(row, model, path, 'simple', checker)
            if claim['checker']['objective'] != old['checked_objective']:
                raise ValueError('legacy panel and replay disagree on the claimed witness')
            claim.update(claimed_status=old['status'], reported_dual_bound=str(old['best_bound']))
            expected.add('reported_bound_excludes_checked_witness')
        facts = assess_evidence(claim, model, relaxation=False, witnesses=witness_index([better]))
        if set(facts['witness_conflicts']) != expected:
            raise ValueError(f'{name} failed historical contradiction control: {facts}')
        controls.append({'instance': name, **facts})
    return controls


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', type=Path, required=True)
    p.add_argument('--verification', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--witness-runs', action='append', default=[], type=Path)
    args = p.parse_args()
    build, out = args.build.resolve(), args.out.resolve()
    receipt = args.verification.resolve()
    candidate = verified(build, receipt)['identity']
    if candidate['project_sanitizers'] or candidate['build_type'] != 'Release':
        raise ValueError('fleet release check needs the verified Release solver')
    config = strict_json((ROOT / 'config/fleet-release.json').read_text())
    suite = ROOT / config['suite']
    if sha(suite) != config['suite_sha256']:
        raise ValueError('fleet suite differs from its lock')
    with suite.open() as stream:
        models = {r['id']:r for r in csv.DictReader(stream, delimiter='\t')}
    if len(models) != 6 or any(sha(ROOT / m['path']) != m['sha256'] for m in models.values()):
        raise ValueError('fleet membership or model contents changed')
    ledgers = [ROOT / path for path in config['witness_ledgers']]
    for path in ledgers:
        if sha(path) != config['witness_ledgers'][str(path.relative_to(ROOT))]:
            raise ValueError(f'archived witness ledger changed: {path}')
    ledgers.extend(path.resolve() for path in args.witness_runs)
    archived, ledger_ids = read_witness_ledgers(ledgers)
    out.mkdir(parents=True, exist_ok=False)
    os.chdir(ROOT)
    frozen = out / 'bin'
    frozen.mkdir()
    for name in TOOLS:
        if sha(shutil.copy2(build / name, frozen / name)) != candidate['binaries'][name]:
            raise ValueError(f'candidate executable changed: {name}')
    checker = frozen / 'milpcheck'
    write(out / 'controls.json', replay_controls(config, models, checker, archived))
    write(out / 'campaign.json', {'identity': candidate, 'verification_sha256': sha(receipt),
          'suite_sha256': sha(suite), 'witness_ledgers': ledger_ids, 'seed': 1, 'time_limit_s': 60})
    solver_file = out / 'solvers.tsv'
    solver_file.write_text('solver\tbinary\tsha256\tversion\n'
                           f"lab\t{frozen/'milpsolve'}\t{candidate['binaries']['milpsolve']}\tverified\n")
    runs = out / 'runs'
    env = {k:v for k,v in os.environ.items() if not k.startswith('LAB_')}
    subprocess.run([frozen/'milpbench','--mip','--suite',suite,'--solvers',solver_file,'--checker',checker,
                    '--out',runs,'--seed','1','--limit','60'], env=env, check=True)
    extra = [arg for path in ledgers for arg in ('--witness-runs',path)]
    subprocess.run([sys.executable, ROOT/'scripts/revalidate.py', runs, '--checker', checker, *extra], check=True)
    subprocess.run([sys.executable, ROOT/'scripts/report.py', runs, '--runs','revalidated.jsonl', *extra], check=True)
    records, _ = read_witness_ledgers([runs/'revalidated.jsonl'])
    if len(records) != len(models):
        raise ValueError('fleet campaign contains missing or duplicate attempts')
    failures = []
    snapshots = []
    for name, model in models.items():
        r = one(records, instance=name)
        facts = r['evaluation']
        status = r['claimed_status']
        snapshots.append({'instance': name, 'model_sha256': r['model_sha256'],
                          **{key: facts[key] for key in (
                              'claimed_status', 'search_route', 'reported_bound', 'checked_objective',
                              'checker_outcome', 'witness_conflicts', 'reference_contradiction',
                              'independently_certified_optimality', 'independently_certified_infeasibility')}})
        if (r['model_sha256'] != model['sha256'] or r['solver_sha256'] != candidate['binaries']['milpsolve']
                or r['runner_sha256'] != candidate['binaries']['milpbench'] or r['checker_sha256'] != sha(checker)
                or r['signal'] or r['watchdog_timeout'] or r['status_conflict'] or not facts['search_route']
                or r['exit_code'] != (0 if status == 'optimal' else 4)
                or status not in ('optimal','time_limit','node_limit','iteration_limit','gap_limit')
                or facts['checker_outcome'] not in ('valid','absent')
                or ('solution.sol' in r['evidence_sha256'] and not facts['checked_feasible'])
                or (status == 'optimal' and not facts['checked_feasible'])
                or (facts['checked_feasible'] and checker_policy_errors(r['checker']))
                or facts['reference_contradiction'] or facts['witness_conflicts']):
            failures.append(name)
    write(out/'release-check.json', {'attempts':len(records), 'failures':failures, 'passed':not failures,
                                    'snapshots':snapshots})
    return int(bool(failures))


if __name__ == '__main__':
    raise SystemExit(main())
