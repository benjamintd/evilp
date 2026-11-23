"""Interpret independent checker evidence for timed runs and revalidations."""
import hashlib
import json
import math
import re
from fractions import Fraction as F
from pathlib import Path

EVALUATION_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

REFERENCE_ABSOLUTE = F(1, 10000)
REFERENCE_RELATIVE = F(1, 10**9)
CHECKER_LIMITS = {'linear_tolerance': F(1, 10000),
                  'integer_tolerance': F(1, 100000),
                  'objective_tolerance': F(1, 10000)}
UNRESOLVED = {'gap_limit', 'unsupported', 'numerical_failure', 'time_limit', 'node_limit',
              'iteration_limit', 'watchdog_timeout', 'error', 'unknown'}


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


def checker_outcome(run):
    status = (run.get('checker') or {}).get('status')
    if status != 'unavailable':
        return status
    code = run.get('checker_exit_code')
    if code == 4:
        return 'unavailable'
    # Older evaluators used 'unavailable' both when no artifacts existed
    # and when the checker timed out or returned an unexpected exit code.
    if code is not None or run.get('checker_argv') or run.get('checker_wall_s', 0) > 0:
        return 'error'
    return 'absent'


def checker_policy_errors(checker):
    reasons = []
    for key, maximum in CHECKER_LIMITS.items():
        try:
            tolerance = rational(checker.get(key))
            if tolerance < 0 or tolerance > maximum:
                reasons.append('checker_tolerance_exceeds_policy_' + key)
        except (ValueError, ZeroDivisionError):
            reasons.append('checker_policy_missing_' + key)
    for key in ['row_violations', 'bound_violations', 'integrality_violations']:
        if type(checker.get(key)) is not int or checker[key] != 0:
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
    return reasons


def read_witness_ledgers(paths):
    rows, identities = [], []
    for path in map(Path, paths):
        content = path.read_bytes()
        rows.extend(strict_json(line) for line in content.splitlines() if line.strip())
        identities.append({'path': str(path.resolve()), 'sha256': hashlib.sha256(content).hexdigest()})
    return rows, identities


def witness_index(rows):
    """Best policy-checked original-MIP incumbent per immutable model and sense."""
    best = {}
    for row in rows:
        model = row.get('model_sha256')
        checker = row.get('checker') or {}
        sense = checker.get('objective_sense')
        if (not isinstance(model, str) or not re.fullmatch('[a-f0-9]{64}', model)
                or row.get('mode') != 'mip' or row.get('relaxation') is not False
                or sense not in ('min', 'max') or checker.get('relaxation') is not False
                or row.get('checker_exit_code') != 0 or checker.get('status') != 'valid'
                or checker.get('checked_feasible') is not True or checker_policy_errors(checker)):
            continue
        try:
            value = rational(checker.get('objective'))
        except (ValueError, ZeroDivisionError):
            continue
        key = (model, sense)
        if key in best and (value >= best[key][0] if sense == 'min' else value <= best[key][0]):
            continue
        source = {field: row.get(field) for field in
                  ('model_sha256', 'instance', 'solver', 'solver_sha256', 'campaign_id', 'seed',
                   'checker_sha256', 'evidence_sha256', 'solution_path')}
        source.update(objective=str(value), objective_sense=sense,
                      tolerances={key: checker[key] for key in CHECKER_LIMITS})
        best[key] = (value, source)
    return best


def excludes(value, incumbent, sense):
    delta = value - incumbent if sense == 'min' else incumbent - value
    return delta > REFERENCE_ABSOLUTE + REFERENCE_RELATIVE * abs(incumbent)


def witness_conflicts(run, objective, index):
    if not index or run.get('mode') != 'mip' or run.get('relaxation') is not False:
        return [], None
    model = run.get('model_sha256')
    sense = (run.get('checker') or {}).get('objective_sense')
    if sense not in ('min', 'max'):
        senses = [s for s in ('min', 'max') if (model, s) in index]
        if len(senses) != 1:
            return [], None
        sense = senses[0]
    if (model, sense) not in index:
        return [], None
    incumbent, source = index[model, sense]
    conflicts = []
    if run.get('claimed_status') == 'infeasible':
        conflicts.append('infeasible_despite_checked_witness')
    if run.get('claimed_status') == 'optimal' and objective is not None and excludes(objective, incumbent, sense):
        conflicts.append('optimal_objective_worse_than_checked_witness')
    if run.get('reported_dual_bound') is not None:
        try:
            bound = rational(run['reported_dual_bound'])
        except (ValueError, ZeroDivisionError):
            bound = None
        if bound is not None and excludes(bound, incumbent, sense):
            conflicts.append('reported_bound_excludes_checked_witness')
    return conflicts, source if conflicts else None


def assess_evidence(run, reference, *, relaxation, witnesses=None):
    """Validate evidence scope before comparing with the requested reference."""
    checker = run.get('checker') or {}
    status = run.get('claimed_status')
    mode = 'relax' if relaxation else 'mip'
    scoped = (run.get('mode', mode) == mode and run.get('relaxation', relaxation) is relaxation
              and checker.get('relaxation') is relaxation)
    if 'sha256' in reference:
        scoped = scoped and run.get('model_sha256') == reference['sha256']
    outcome = checker_outcome(run)
    if outcome == 'valid' and (not scoped or run.get('checker_exit_code') != 0):
        outcome = 'error'
    valid = outcome == 'valid'
    # The checker can reject the proposed primal while independently
    # validating a separately supplied dual bound on the same model.
    bound_checked = valid or (scoped and outcome == 'invalid' and run.get('checker_exit_code') == 1)

    def number(value):
        try:
            return rational(value)
        except (ValueError, ZeroDivisionError):
            return None

    objective = number(checker.get('objective')) if valid and checker.get('checked_feasible') is True else None
    bound = number(checker.get('certified_bound')) if bound_checked and checker.get('bound_status') == 'valid' else None
    expected = reference.get('reference_status')
    target = number(reference.get('reference_objective')) if expected in ('=opt=', '=best=') else None
    feasible = objective is not None
    match = (feasible and target is not None and
             abs(objective - target) <= REFERENCE_ABSOLUTE + REFERENCE_RELATIVE * abs(target))
    better = (feasible and target is not None and not match and expected == '=opt=' and
              ((checker.get('objective_sense') == 'min' and objective < target) or
               (checker.get('objective_sense') == 'max' and objective > target)))
    infeasibility_proof = (valid and checker.get('certificate_status') == 'valid'
                           and checker.get('certifies_infeasibility') is True)
    contradiction = ((infeasibility_proof and (feasible or expected in ('=opt=', '=best='))) or
                     (status == 'infeasible' and (feasible or expected in ('=opt=', '=best='))) or
                     (feasible and expected == '=inf=') or
                     (status == 'optimal' and feasible and expected == '=opt=' and target is not None and not match)
                     or better)
    sense = checker.get('objective_sense')
    if target is not None and sense in ('min', 'max'):
        raw_bound = number(run.get('reported_dual_bound'))
        if ((raw_bound is not None and excludes(raw_bound, target, sense)) or
                (status == 'optimal' and feasible and excludes(objective, target, sense))):
            contradiction = True
    conflicts, source = ([], None) if relaxation else witness_conflicts(run, objective, witnesses)
    gap = None
    raw_bound = number(run.get('reported_dual_bound'))
    if (not relaxation and feasible and raw_bound is not None and not conflicts and not contradiction
            and run.get('signal') == 0 and run.get('watchdog_timeout') is False
            and run.get('status_conflict') is False and run.get('exit_code') in ((0,) if status == 'optimal' else (0, 4))
            and run.get('bound_origin') not in (None, 'none')
            and run.get('bound_scope') == 'original_mip' and sense in ('min', 'max')
            and status in ('optimal', 'time_limit', 'node_limit', 'iteration_limit', 'gap_limit')
            and not checker_policy_errors(checker)):
        difference = objective - raw_bound if sense == 'min' else raw_bound - objective
        if difference == 0:
            gap = F(0)
        elif difference > 0 and objective * raw_bound > 0:
            gap = difference / min(abs(objective), abs(raw_bound))
    return {'claimed_status': status, 'reported_bound': run.get('reported_dual_bound'),
            **{key: run.get(key) for key in ('search_route', 'closure_reason', 'bound_origin', 'bound_scope')},
            'reported_relative_gap': str(gap) if gap is not None else None,
            'evidence_scope': mode if scoped else None, 'checker_outcome': outcome,
            'checked_feasible': feasible, 'witness_conflicts': conflicts, 'contradicting_witness': source,
            'checked_objective': str(objective) if feasible else None,
            'reference_objective_agreement': match, 'reference_contradiction': contradiction,
            'certified_bound': str(bound) if bound is not None else None,
            'independently_certified_optimality': (feasible and bound == objective and
                                                  checker.get('certifies_optimality') is True),
            'independently_certified_infeasibility': infeasibility_proof}


def classify(model, run, witnesses=None):
    reasons = []
    status = run.get('claimed_status')
    checker = run.get('checker') or {}
    facts = assess_evidence(run, model, relaxation=False, witnesses=witnesses)
    reasons.extend(facts['witness_conflicts'])
    if facts['reference_contradiction']:
        reasons.append('reference_contradiction')
    feasible = facts['checked_feasible']
    if checker.get('checked_feasible') is True and checker.get('relaxation') is not False:
        reasons.append('original_MILP_checker_scope_required')
    independent_infeasible = facts['independently_certified_infeasibility']
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
    if facts['checker_outcome'] in ('invalid', 'error', 'unsupported'):
        reasons.append('checker_' + facts['checker_outcome'])
    # Some adapters emit a header-only solution file for an infeasible model.
    # The checker explicitly recognizes this format and returns no_solution (4).
    # An absent primal is expected only under the declared status-only contract.
    empty_primal = (model['reference_status'] == '=inf=' and status == 'infeasible' and
                    model.get('validation_contract') == 'reference_status_agreement' and
                    run.get('checker_exit_code') == 4 and checker.get('status') == 'no_solution' and
                    checker.get('checked_feasible') is False and checker.get('objective') is None and
                    checker.get('certificate_status') == 'absent' and
                    checker.get('certifies_infeasibility') is False and
                    checker.get('certifies_optimality') is False)
    if run.get('checker_exit_code') not in (None, 0) and not empty_primal:
        reasons.append('checker_process_error')
    if run.get('reference_status') != model['reference_status']:
        reasons.append('reference_status_mismatch')
    expected = model['reference_status']
    match = facts['reference_objective_agreement']
    if expected in ('=opt=', '=best='):
        try:
            reference = rational(model['reference_objective'])
            if rational(run.get('reference_objective')) != reference:
                reasons.append('reference_objective_mismatch')
        except (ValueError, ZeroDivisionError):
            reasons.append('malformed_reference_or_checked_objective')
    if independent_infeasible and (feasible or expected in ('=opt=', '=best=')):
        reasons.append('infeasibility_proof_conflicts_with_feasible_evidence_or_reference')
    if expected == '=opt=':
        if status != 'optimal':
            reasons.append('optimal_completion_required')
        if not feasible or run.get('checker_exit_code') != 0:
            reasons.append('independently_checked_feasible_witness_required')
        if not match:
            reasons.append('locked_reference_not_matched')
        if feasible:
            reasons.extend(checker_policy_errors(checker))
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
    return {**facts, 'completed_correctly': not reasons, 'reasons': sorted(set(reasons)),
            'reference_status_agreement_is_not_an_independent_proof': True}
