#!/usr/bin/env python3
"""Full MIPLIB completion and paired HiGHS speed gate, with no optimization.

The locked manifest must contain exactly 240 distinct models. Optimal status,
reference agreement, an independently checked feasible point and a proof remain
separate fields. Unknown references cannot pass. Known infeasible cases require
their explicit validation_contract or an original-model independent proof.

Speed evidence must declare a serial, quiescent paired campaign and bind every
record to its campaign, hardware, budget and engine configuration. Historical
or concurrent-engineering inventories cannot open the speed gate. Wins require
strictly smaller solver process wall time and valid completion by both engines;
ties, missing results and limits never win. The denominator is always 240.
"""
import argparse
import collections
import csv
import hashlib
import json
import pathlib
import re
import sys

from run_evaluation import (EVALUATION_SHA256, REFERENCE_ABSOLUTE, REFERENCE_RELATIVE, UNRESOLVED,
                            classify, rational, strict_json, read_witness_ledgers, witness_index)
CORPUS_SIZE = 240
REQUIRED_WINS = 120
SHA256 = re.compile(r'[0-9a-f]{64}\Z')
MEASUREMENT = 'solver_process_wall_excluding_checker'
BUDGET_KEYS = ('time_limit_s', 'watchdog_grace_s', 'seed', 'requested_threads',
               'memory_limit_bytes', 'node_limit', 'iteration_limit')


def digest(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def validate_manifest(manifest):
    errors = []
    instances = manifest.get('instances')
    if not isinstance(instances, list):
        return {}, ['manifest.instances must be a complete list']
    if manifest.get('count') != CORPUS_SIZE or len(instances) != CORPUS_SIZE:
        errors.append('manifest must declare and contain all 240 instances')
    models = {}
    for index, model in enumerate(instances):
        if not isinstance(model, dict):
            errors.append(f'manifest instance {index} is not an object')
            continue
        name = model.get('id')
        if not isinstance(name, str) or not name or name in models:
            errors.append(f'missing or duplicate manifest instance: {name!r}')
            continue
        models[name] = model
        if not isinstance(model.get('sha256'), str) or not SHA256.fullmatch(model['sha256']):
            errors.append(name + ': valid locked model SHA256 required')
        if not isinstance(model.get('reference_status'), str):
            errors.append(name + ': locked reference status required')
        if model.get('reference_status') in ('=opt=', '=best='):
            try:
                rational(model.get('reference_objective'))
            except (ValueError, ZeroDivisionError):
                errors.append(name + ': finite locked reference objective required')
    if len(models) != CORPUS_SIZE:
        errors.append('manifest membership is not exactly 240 unique instances')
    if len({model.get('sha256') for model in models.values()}) != len(models):
        errors.append('duplicate model hashes cannot replace distinct corpus members')
    return models, errors


def index_runs(runs, models, solver, expected_engine=None):
    indexed, errors = {}, []
    identities = {}
    for run in runs:
        if not isinstance(run, dict):
            errors.append(solver + ': run must be a JSON object')
            continue
        name = run.get('instance')
        if run.get('solver') != solver:
            errors.append(solver + ': unexpected solver label ' + str(run.get('solver')))
        if name not in models:
            errors.append(solver + ': unexpected model ' + str(name))
        if name in indexed:
            errors.append(solver + ': duplicate attempt for ' + str(name))
        for key, description in (('solver_sha256', 'solver binaries'),
                                 ('configuration_sha256', 'solver configurations'),
                                 ('runner_sha256', 'runner binaries'), ('checker_sha256', 'checker binaries')):
            observed = run.get(key)
            if not isinstance(observed, str) or not SHA256.fullmatch(observed):
                errors.append(f'{solver}/{name}: missing or invalid {key}')
            elif observed != identities.setdefault(key, observed):
                errors.append(f'{solver}: attempts use different {description}')
            if expected_engine is not None and (key in expected_engine or key in ('solver_sha256', 'configuration_sha256')):
                if expected_engine.get(key) != observed:
                    errors.append(f'{solver}/{name}: differs from declared {key}')
        indexed[name] = run
    missing = sorted(set(models)-set(indexed))
    if missing:
        errors.append(solver + ': missing ' + str(len(missing)) + ' model attempts')
    if len(runs) != len(models) or len(indexed) != len(models):
        errors.append(f'{solver}: exactly {len(models)} unique attempts required')
    return indexed, sorted(set(errors)), missing


def timing_contract(campaign, candidate, baseline, candidate_solver, baseline_solver):
    errors = []
    if not isinstance(campaign, dict):
        return ['paired timing campaign metadata is missing']
    for key in ['fair_timing_campaign', 'serial_solver_processes', 'quiescent']:
        if campaign.get(key) is not True:
            errors.append(key + ' must be explicitly true')
    if campaign.get('coverage_exploration') is not False:
        errors.append('coverage exploration cannot establish speed ranking')
    if not campaign.get('campaign_id') or not campaign.get('hardware_id'):
        errors.append('paired campaign and hardware identities required')
    if campaign.get('measurement') != MEASUREMENT:
        errors.append('matched solver-process wall measurement required')
    budget = campaign.get('budget')
    if not isinstance(budget, dict) or any(key not in budget for key in BUDGET_KEYS):
        errors.append('complete common budget/caps declaration required')
        budget = {}
    else:
        if type(budget['requested_threads']) is not int or budget['requested_threads'] != 1:
            errors.append('single-thread matched campaign required')
        for key in ['time_limit_s', 'watchdog_grace_s']:
            try:
                if rational(budget[key]) <= 0:
                    errors.append('positive finite ' + key + ' required')
            except (ValueError, ZeroDivisionError):
                errors.append('positive finite ' + key + ' required')
        if not isinstance(budget['seed'], int) or isinstance(budget['seed'], bool) or budget['seed'] <= 0:
            errors.append('positive integer seed required')
        for key in ['memory_limit_bytes', 'node_limit', 'iteration_limit']:
            if budget[key] is not None and (type(budget[key]) is not int or budget[key] <= 0):
                errors.append('positive integer or explicit unlimited ' + key + ' required')
    engines = campaign.get('engines') or {}
    for name in set(candidate) & set(baseline):
        if candidate[name].get('model_sha256') != baseline[name].get('model_sha256'):
            errors.append(str(name) + ': candidate and HiGHS inputs differ')
    for label, records in [(candidate_solver, candidate), (baseline_solver, baseline)]:
        engine = engines.get(label) or {}
        for key in ['solver_sha256', 'configuration_sha256']:
            if not isinstance(engine.get(key), str) or not SHA256.fullmatch(engine[key]):
                errors.append(label + ': frozen ' + key + ' required')
        for name, run in records.items():
            for key in ('runner_sha256', 'checker_sha256'):
                expected = campaign.get(key)
                if not isinstance(expected, str) or not SHA256.fullmatch(expected) or run.get(key) != expected:
                    errors.append(f'{label}/{name}: unmatched {key}')
            if run.get('sleep_affected') is not False:
                errors.append(label + '/' + str(name) + ': missing or failed host-suspension check')
            if (run.get('campaign_id') != campaign.get('campaign_id') or
                    run.get('hardware_id') != campaign.get('hardware_id') or
                    run.get('measurement') != MEASUREMENT or run.get('budget') != budget or
                    run.get('solver_sha256') != engine.get('solver_sha256') or
                    run.get('configuration_sha256') != engine.get('configuration_sha256')):
                errors.append(label + '/' + str(name) + ': unmatched campaign/hardware/budget/configuration')
            for key in ['time_limit_s', 'watchdog_grace_s', 'seed', 'requested_threads']:
                if run.get(key) != budget.get(key):
                    errors.append(label + '/' + str(name) + ': run differs from declared ' + key)
    return sorted(set(errors))


def evaluate(manifest, candidate_runs, baseline_runs=None, campaign=None,
             candidate_solver='lab', baseline_solver='highs', panel=None, candidate_engine=None, witness_runs=()):
    models, manifest_errors = validate_manifest(manifest)
    locked_count = len(models)
    if panel is not None:
        selected = {row['id']: row for row in panel}
        if not selected or len(selected) != len(panel) or not selected.keys() <= models.keys():
            raise ValueError('panel must contain distinct members of the locked corpus')
        for name, row in selected.items():
            for key in ('sha256', 'reference_status', 'reference_objective'):
                if row.get(key) != models[name].get(key):
                    raise ValueError(f'{name}: panel differs from locked {key}')
        models = {name: models[name] for name in selected}
    engines = (campaign or {}).get('engines', {})
    declared = candidate_engine or engines.get(candidate_solver)
    if declared is not None and campaign:
        declared = {**declared, **{key: campaign.get(key) for key in ('runner_sha256', 'checker_sha256')}}
    candidate, candidate_errors, missing = index_runs(candidate_runs, models, candidate_solver, declared)
    if declared is None:
        candidate_errors.append('completion requires a frozen candidate configuration or campaign')
    witnesses = witness_index([*candidate_runs, *(baseline_runs or ()), *witness_runs])
    classifications = {name: classify(model, candidate[name], witnesses) for name, model in models.items() if name in candidate}
    complete = sum(r['completed_correctly'] for r in classifications.values())
    completion_passed = panel is None and not manifest_errors and not candidate_errors and complete == CORPUS_SIZE
    timing_errors, baseline_errors, comparisons = [], [], []
    baseline, baseline_classifications = {}, {}
    if baseline_runs is None:
        timing_errors.append('paired HiGHS attempts are missing; historical baselines are not substituted')
    else:
        baseline, baseline_errors, _ = index_runs(baseline_runs, models, baseline_solver, engines.get(baseline_solver))
        baseline_classifications = {name: classify(model, baseline[name], witnesses) for name, model in models.items() if name in baseline}
        timing_errors = timing_contract(campaign, candidate, baseline, candidate_solver, baseline_solver)
    comparable = not manifest_errors and not candidate_errors and not baseline_errors and not timing_errors
    wins = ties = losses = eligible = valid_solved_pairs = 0
    for name in models:
        valid_pair = (name in classifications and name in baseline_classifications and
                      classifications[name]['completed_correctly'] and baseline_classifications[name]['completed_correctly'])
        valid_solved_pairs += bool(valid_pair)
        comparison = {'instance': name, 'validly_completed_by_both': bool(valid_pair),
                      'timing_eligible': False, 'candidate_strictly_faster': False}
        if comparable and valid_pair:
            try:
                a, b = rational(candidate[name].get('wall_s')), rational(baseline[name].get('wall_s'))
                if a <= 0 or b <= 0:
                    raise ValueError('positive duration required')
                eligible += 1
                wins += a < b
                ties += a == b
                losses += a > b
                comparison.update(timing_eligible=True, candidate_strictly_faster=a < b,
                                  candidate_wall_s=str(a), highs_wall_s=str(b))
            except (ValueError, ZeroDivisionError):
                comparison['reason'] = 'positive finite wall durations required'
        comparisons.append(comparison)
    blockers = collections.Counter(reason for result in classifications.values() for reason in result['reasons'])
    unresolved_references = sorted(name for name, model in models.items()
                                   if model.get('reference_status') not in ('=opt=', '=inf='))
    speed_passed = completion_passed and comparable and wins >= REQUIRED_WINS
    return {'schema': 1, 'kind': 'Development panel' if panel is not None else 'Full 240-instance completion and matched paired-speed gate',
            'scope': ('Selected locked models only. A panel cannot open the full gate.' if panel is not None else
                      'Full MIPLIB corpus. Other regression suites remain separate prerequisites.'),
            'manifest_count': locked_count, 'selected_count': len(models),
            'panel_records_complete': (panel is not None and baseline_runs is not None and
                                       not (manifest_errors or candidate_errors or baseline_errors)), 'required_instances': CORPUS_SIZE,
            'manifest_errors': manifest_errors, 'candidate_attempt_errors': candidate_errors,
            'missing_candidate_instances': missing, 'completed_correctly': complete,
            'completion_gate_passed': completion_passed, 'speed_gate_passed': speed_passed,
            'passed': completion_passed and speed_passed,
            'reference_policy': {'absolute': str(REFERENCE_ABSOLUTE), 'relative': str(REFERENCE_RELATIVE),
                                 'agreement_is_not_independent_proof': True},
            'reference_status_counts': dict(collections.Counter(m.get('reference_status') for m in models.values())),
            'unresolved_reference_instances': unresolved_references,
            'candidate_status_counts': dict(collections.Counter(r.get('claimed_status') for r in candidate_runs)),
            'candidate_blocker_counts': dict(blockers.most_common()), 'candidate_results': classifications,
            'baseline_results': baseline_classifications, 'baseline_attempt_errors': baseline_errors, 'timing_errors': timing_errors,
            'timing_comparable': comparable, 'paired_valid_solved_count': valid_solved_pairs,
            'timing_eligible_pairs': eligible,
            'strict_candidate_wins': wins, 'ties': ties, 'losses': losses,
            'speed_denominator': CORPUS_SIZE, 'required_strict_wins': REQUIRED_WINS,
            'speed_comparisons': comparisons,
            'never_counts_as_pass': sorted(UNRESOLVED),
            'proof_scope': 'Individual LP-box proofs do not establish a complete MILP tree proof.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=pathlib.Path, default=pathlib.Path('config/miplib.lock.json'))
    parser.add_argument('--runs', type=pathlib.Path, required=True)
    parser.add_argument('--baseline-runs', type=pathlib.Path)
    parser.add_argument('--witness-runs', action='append', type=pathlib.Path, default=[],
                        help='include checked incumbents from another run ledger (repeatable)')
    parser.add_argument('--campaign', type=pathlib.Path)
    parser.add_argument('--candidate-config', type=pathlib.Path, help='runner-produced candidate config for standalone completion')
    parser.add_argument('--panel-suite', type=pathlib.Path, help='score a declared subset without opening the full gate')
    parser.add_argument('--candidate', default='lab')
    parser.add_argument('--baseline', default='highs')
    parser.add_argument('--out', type=pathlib.Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('--out must be a fresh file')
    try:
        manifest = strict_json(args.manifest.read_text())
        runs = [strict_json(line) for line in args.runs.read_text().splitlines() if line.strip()]
        candidate = [run for run in runs if run.get('solver') == args.candidate]
        evidence = [r for r in runs if r.get('solver') != args.candidate]
        baseline = None
        if args.baseline_runs:
            baseline = [strict_json(line) for line in args.baseline_runs.read_text().splitlines() if line.strip()]
            evidence.extend(r for r in baseline if r.get('solver') != args.baseline)
            baseline = [run for run in baseline if run.get('solver') == args.baseline]
        campaign = strict_json(args.campaign.read_text()) if args.campaign else None
        engine = panel = None
        if args.candidate_config:
            config = strict_json(args.candidate_config.read_text())
            engine = {'solver_sha256': config['solver_sha256'], 'configuration_sha256': digest(args.candidate_config)}
            if campaign and engine != campaign['engines'].get(args.candidate):
                raise ValueError('candidate configuration differs from campaign')
            if (config['adapter'] != args.candidate or config['mode'] != 'mip' or
                    (campaign and config['runner_sha256'] != campaign.get('runner_sha256'))):
                raise ValueError('candidate configuration has an incompatible runner, mode or adapter')
            engine['runner_sha256'] = config['runner_sha256']
        if args.panel_suite:
            if not campaign or campaign.get('suite_sha256') != digest(args.panel_suite):
                raise ValueError('panel suite must match the frozen campaign suite')
            with args.panel_suite.open() as stream:
                panel = list(csv.DictReader(stream, delimiter='\t'))
        archived, ledgers = read_witness_ledgers(args.witness_runs)
        evidence.extend(archived)
        result = evaluate(manifest, candidate, baseline, campaign, args.candidate, args.baseline, panel, engine, evidence)
        result['witness_ledgers'] = ledgers
        result['inputs'] = {name: {'path': str(path.resolve()), 'sha256': digest(path)}
                            for name, path in [('manifest', args.manifest), ('runs', args.runs),
                                               ('baseline_runs', args.baseline_runs), ('campaign', args.campaign),
                                               ('candidate_config', args.candidate_config), ('panel_suite', args.panel_suite)] if path}
        result['driver_sha256'] = digest(__file__)
        result['evaluation_sha256'] = EVALUATION_SHA256
    except (ValueError, OSError, KeyError, TypeError, AttributeError) as error:
        result = {'schema': 1, 'passed': False, 'input_error': str(error)}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: result.get(key) for key in ['passed', 'completion_gate_passed',
          'completed_correctly', 'strict_candidate_wins', 'speed_gate_passed', 'input_error']}))
    if 'input_error' in result:
        return 2
    return 0 if result['passed'] or result['panel_records_complete'] else 1


if __name__ == '__main__':
    sys.exit(main())
