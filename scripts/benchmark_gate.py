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
import fractions
import hashlib
import json
import math
import pathlib
import re
import sys

F = fractions.Fraction
CORPUS_SIZE = 240
REQUIRED_WINS = 120
SHA256 = re.compile(r'[0-9a-f]{64}\Z')
REFERENCE_ABSOLUTE = F(1, 10000)
REFERENCE_RELATIVE = F(1, 10**9)
CHECKER_LIMITS = {'linear_tolerance': F(1, 10000),
                  'integer_tolerance': F(1, 100000),
                  'objective_tolerance': F(1, 10000)}
UNRESOLVED = {'unsupported', 'numerical_failure', 'time_limit', 'node_limit',
              'iteration_limit', 'watchdog_timeout', 'error', 'unknown'}
MEASUREMENT = 'solver_process_wall_excluding_checker'
BUDGET_KEYS = ('time_limit_s', 'watchdog_grace_s', 'seed', 'requested_threads',
               'memory_limit_bytes', 'node_limit', 'iteration_limit')


def digest(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def unique_object(pairs):
    out = {}
    for key, value in pairs:
        if key in out:
            raise ValueError('duplicate JSON key: ' + key)
        out[key] = value
    return out


def strict_json(text):
    def invalid(value):
        raise ValueError('nonfinite JSON number: ' + value)
    def finite_float(value):
        parsed = float(value)
        if not math.isfinite(parsed):
            invalid(value)
        return parsed
    return json.loads(text, object_pairs_hook=unique_object, parse_constant=invalid,
                      parse_float=finite_float)


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        raise ValueError('finite rational value required')
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError('nonfinite value')
    return F(str(value))


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


def classify(model, run):
    reasons = []
    status = run.get('claimed_status')
    checker = run.get('checker') or {}
    feasible = checker.get('status') == 'valid' and checker.get('checked_feasible') is True
    independent_optimum = checker.get('certifies_optimality') is True
    independent_infeasible = (checker.get('certificate_status') == 'valid' and
                              checker.get('certifies_infeasibility') is True)
    if run.get('model_sha256') != model['sha256']:
        reasons.append('model_hash_mismatch')
    if run.get('mode') != 'mip' or run.get('relaxation') is not False:
        reasons.append('MILP_mode_required')
    if run.get('exit_code') != 0 or run.get('signal') != 0 or run.get('watchdog_timeout') is not False:
        reasons.append('process_did_not_complete_cleanly')
    if run.get('status_conflict'):
        reasons.append('conflicting_solver_status')
    if status in UNRESOLVED:
        reasons.append('unresolved_' + str(status))
    if checker.get('status') in ('invalid', 'error', 'unsupported'):
        reasons.append('checker_' + checker['status'])
    if run.get('checker_exit_code') not in (None, 0):
        reasons.append('checker_process_error')
    if run.get('reference_status') != model['reference_status']:
        reasons.append('reference_status_mismatch')
    expected = model['reference_status']
    match, objective = False, None
    if expected in ('=opt=', '=best='):
        try:
            reference = rational(model['reference_objective'])
            if rational(run.get('reference_objective')) != reference:
                reasons.append('reference_objective_mismatch')
            if feasible:
                objective = rational(checker.get('objective'))
                match = abs(objective-reference) <= REFERENCE_ABSOLUTE + REFERENCE_RELATIVE*abs(reference)
        except (ValueError, ZeroDivisionError):
            reasons.append('malformed_reference_or_checked_objective')
    if expected == '=opt=':
        if status != 'optimal':
            reasons.append('optimal_completion_required')
        if not feasible or run.get('checker_exit_code') != 0:
            reasons.append('independently_checked_feasible_witness_required')
        if not match:
            reasons.append('locked_reference_not_matched')
        if feasible:
            for key, maximum in CHECKER_LIMITS.items():
                try:
                    tolerance = rational(checker.get(key))
                    if tolerance < 0 or tolerance > maximum:
                        reasons.append('checker_tolerance_exceeds_policy_' + key)
                except (ValueError, ZeroDivisionError):
                    reasons.append('checker_policy_missing_' + key)
            for key in ['row_violations', 'bound_violations', 'integrality_violations']:
                if checker.get(key) != 0:
                    reasons.append('checker_violation_or_missing_count_' + key)
            for key, tolerance in [('max_row_violation', 'linear_tolerance'),
                                   ('max_bound_violation', 'linear_tolerance'),
                                   ('max_integrality_violation', 'integer_tolerance'),
                                   ('objective_error', 'objective_tolerance')]:
                try:
                    value = rational(checker.get(key))
                    if value < 0 or value > rational(checker.get(tolerance)):
                        reasons.append('checker_report_inconsistent_' + key)
                except (ValueError, ZeroDivisionError):
                    reasons.append('checker_measurement_missing_' + key)
    elif expected == '=inf=':
        if status != 'infeasible':
            reasons.append('infeasible_completion_required')
        contract = model.get('validation_contract', 'independent_original_model_infeasibility')
        if contract == 'reference_status_agreement':
            # Explicitly declared reference-based validation is not a proof.
            if not model.get('validation_basis'):
                reasons.append('declared_infeasibility_validation_basis_required')
        elif contract == 'independent_original_model_infeasibility':
            if not independent_infeasible or run.get('checker_exit_code') != 0:
                reasons.append('independent_original_model_infeasibility_required')
        else:
            reasons.append('unsupported_infeasibility_validation_contract')
        if feasible:
            reasons.append('infeasibility_contradicted_by_checked_witness')
    else:
        reasons.append('unresolved_expected_outcome_' + str(expected))
    return {'completed_correctly': not reasons, 'reasons': sorted(set(reasons)),
            'claimed_status': status, 'checked_feasible': feasible,
            'reference_objective_agreement': match,
            'checked_objective': str(objective) if objective is not None else None,
            'independently_certified_optimality': independent_optimum,
            'independently_certified_infeasibility': independent_infeasible,
            'reference_status_agreement_is_not_an_independent_proof': True}


def index_runs(runs, models, solver):
    indexed, errors = {}, []
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
        indexed[name] = run
    missing = sorted(set(models)-set(indexed))
    if missing:
        errors.append(solver + ': missing ' + str(len(missing)) + ' model attempts')
    if len(runs) != CORPUS_SIZE or len(indexed) != CORPUS_SIZE:
        errors.append(solver + ': exactly 240 unique attempts required')
    return indexed, errors, missing


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
             candidate_solver='lab', baseline_solver='highs'):
    models, manifest_errors = validate_manifest(manifest)
    candidate, candidate_errors, missing = index_runs(candidate_runs, models, candidate_solver)
    classifications = {name: classify(model, candidate[name]) for name, model in models.items() if name in candidate}
    complete = sum(r['completed_correctly'] for r in classifications.values())
    completion_passed = not manifest_errors and not candidate_errors and complete == CORPUS_SIZE
    timing_errors, baseline_errors, comparisons = [], [], []
    baseline, baseline_classifications = {}, {}
    if baseline_runs is None:
        timing_errors.append('paired HiGHS attempts are missing; historical baselines are not substituted')
    else:
        baseline, baseline_errors, _ = index_runs(baseline_runs, models, baseline_solver)
        baseline_classifications = {name: classify(model, baseline[name]) for name, model in models.items() if name in baseline}
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
    return {'schema': 1, 'kind': 'Full 240-instance completion and matched paired-speed gate',
            'scope': 'MIPLIB corpus only; all other existing regression suites remain separate required prerequisites.',
            'manifest_count': len(models), 'required_instances': CORPUS_SIZE,
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
            'baseline_attempt_errors': baseline_errors, 'timing_errors': timing_errors,
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
    parser.add_argument('--campaign', type=pathlib.Path)
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
        baseline = None
        if args.baseline_runs:
            baseline = [strict_json(line) for line in args.baseline_runs.read_text().splitlines() if line.strip()]
            baseline = [run for run in baseline if run.get('solver') == args.baseline]
        campaign = strict_json(args.campaign.read_text()) if args.campaign else None
        result = evaluate(manifest, candidate, baseline, campaign, args.candidate, args.baseline)
        result['inputs'] = {name: {'path': str(path.resolve()), 'sha256': digest(path)}
                            for name, path in [('manifest', args.manifest), ('runs', args.runs),
                                               ('baseline_runs', args.baseline_runs), ('campaign', args.campaign)] if path}
        result['driver_sha256'] = digest(__file__)
    except (ValueError, OSError, KeyError, TypeError, AttributeError) as error:
        result = {'schema': 1, 'passed': False, 'input_error': str(error)}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: result.get(key) for key in ['passed', 'completion_gate_passed',
          'completed_correctly', 'strict_candidate_wins', 'speed_gate_passed', 'input_error']}))
    return 2 if 'input_error' in result else 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
