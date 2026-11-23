#!/usr/bin/env python3
"""Adversarial accounting contracts; no optimization solver is invoked."""
import copy
import json
import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'scripts'))
from benchmark_gate import evaluate, strict_json, MEASUREMENT


def fixture():
    models = [{'id': f'm{i:03d}', 'sha256': f'{i:064x}', 'reference_status': '=opt=',
               'reference_objective': '1'} for i in range(240)]
    manifest = {'count': 240, 'instances': models}
    budget = {'time_limit_s': 10, 'watchdog_grace_s': 15, 'seed': 1, 'requested_threads': 1,
              'memory_limit_bytes': None, 'node_limit': None, 'iteration_limit': None}
    campaign = {'campaign_id': 'paired-new-campaign', 'hardware_id': 'same-host-fingerprint',
                'runner_sha256': '1'*64, 'checker_sha256': '2'*64,
                'fair_timing_campaign': True, 'serial_solver_processes': True, 'quiescent': True,
                'coverage_exploration': False, 'measurement': MEASUREMENT, 'budget': budget,
                'engines': {name: {'solver_sha256': char*64, 'configuration_sha256': config*64}
                            for name, char, config in [('lab', 'a', 'b'), ('highs', 'c', 'd')]}}
    def runs(solver):
        return [{'instance': model['id'], 'solver': solver, 'model_sha256': model['sha256'],
                 'mode': 'mip', 'relaxation': False, 'claimed_status': 'optimal', 'exit_code': 0, 'sleep_affected': False,
                 'runner_sha256': '1'*64, 'checker_sha256': '2'*64,
                 'signal': 0, 'watchdog_timeout': False, 'status_conflict': False,
                 'reference_status': '=opt=', 'reference_objective': '1', 'checker_exit_code': 0,
                 'checker': {'status': 'valid', 'checked_feasible': True, 'relaxation': False, 'objective': '1',
                             'linear_tolerance': '1e-4', 'integer_tolerance': '1e-5',
                             'objective_tolerance': '1e-4', 'certifies_optimality': False,
                             'row_violations': 0, 'bound_violations': 0, 'integrality_violations': 0,
                             'max_row_violation': '0', 'max_bound_violation': '0',
                             'max_integrality_violation': '0', 'objective_error': '0'},
                 'wall_s': 1 if solver == 'lab' else 2, 'campaign_id': campaign['campaign_id'],
                 'hardware_id': campaign['hardware_id'], 'measurement': MEASUREMENT,
                 'budget': copy.deepcopy(budget), **budget, **campaign['engines'][solver]}
                for model in models]
    return manifest, runs('lab'), runs('highs'), campaign


class Gate(unittest.TestCase):
    def setUp(self):
        self.manifest, self.lab, self.highs, self.campaign = fixture()

    def gate(self):
        return evaluate(self.manifest, self.lab, self.highs, self.campaign)

    def test_checked_incumbent_overrules_optimal_and_speed_claims(self):
        for sense, value in [('min', '0.9998'), ('max', '1.0002')]:
            witness = copy.deepcopy(self.highs[0])
            witness.update(claimed_status='time_limit', exit_code=4, seed=42, instance='renamed')
            witness['checker'].update(objective=value, objective_sense=sense)
            self.lab[0]['checker']['objective_sense'] = sense
            self.lab[0]['reported_dual_bound'] = '1'
            result = evaluate(self.manifest, self.lab, self.highs, self.campaign, witness_runs=[witness])
            facts = result['candidate_results']['m000']
            self.assertCountEqual(facts['witness_conflicts'], ['reported_bound_excludes_checked_witness',
                                                              'optimal_objective_worse_than_checked_witness'])
            self.assertEqual(facts['contradicting_witness']['seed'], 42)
            self.assertFalse(facts['independently_certified_optimality'])
            self.assertFalse(result['completion_gate_passed'])
            self.assertEqual(result['strict_candidate_wins'], 239)

    def test_witness_eligibility_and_limits(self):
        from run_evaluation import assess_evidence, witness_index
        row = self.lab[0]
        row['checker']['objective_sense'] = 'min'
        witness = copy.deepcopy(row)
        witness['checker']['objective'] = '0'
        for mutation in [dict(model_sha256='f'*64), dict(mode='relax'), dict(checker_exit_code=1),
                         dict(checker={**witness['checker'], 'relaxation': True}),
                         dict(checker={**witness['checker'], 'linear_tolerance': '0.1'}),
                         dict(checker={**witness['checker'], 'max_row_violation': '2'})]:
            facts = assess_evidence(row, row, relaxation=False, witnesses=witness_index([{**witness, **mutation}]))
            self.assertFalse(facts['witness_conflicts'], mutation)
        row.update(claimed_status='time_limit', reported_dual_bound='0.0001')
        index = witness_index([witness])
        self.assertFalse(assess_evidence(row, row, relaxation=False, witnesses=index)['witness_conflicts'])
        row['reported_dual_bound'] = '0.00010000000000000001'
        self.assertEqual(assess_evidence(row, row, relaxation=False, witnesses=index)['witness_conflicts'],
                         ['reported_bound_excludes_checked_witness'])
        row.update(claimed_status='infeasible', checker={}, reported_dual_bound=None)
        self.assertEqual(assess_evidence(row, row, relaxation=False, witnesses=index)['witness_conflicts'],
                         ['infeasible_despite_checked_witness'])

    def test_all_240_complete_and_120_strict_wins_pass(self):
        for record in self.lab[120:]:
            record['wall_s'] = 3
        result = self.gate()
        self.assertTrue(result['passed'])
        self.assertEqual(result['strict_candidate_wins'], 120)
        self.assertEqual(result['speed_denominator'], 240)
        self.assertFalse(result['candidate_results']['m000']['independently_certified_optimality'])

    def test_mixed_candidate_build_or_configuration_blocks_standalone_completion(self):
        for key in ('solver_sha256', 'configuration_sha256', 'runner_sha256', 'checker_sha256'):
            with self.subTest(key=key):
                manifest, candidate, _, _ = fixture()
                candidate[0][key] = 'f' * 64
                self.assertFalse(evaluate(manifest, candidate, candidate_engine=self.campaign['engines']['lab'])['completion_gate_passed'])
                del candidate[0][key]
                self.assertFalse(evaluate(manifest, candidate, candidate_engine=self.campaign['engines']['lab'])['completion_gate_passed'])

    def test_standalone_completion_requires_declared_build(self):
        self.assertFalse(evaluate(self.manifest, self.lab)['completion_gate_passed'])
        result = evaluate(self.manifest, self.lab, candidate_engine=self.campaign['engines']['lab'])
        self.assertTrue(result['completion_gate_passed'])
        self.assertFalse(result['speed_gate_passed'])

    def test_panel_membership_and_accounting(self):
        panel = self.manifest['instances'][:3]
        self.lab[0]['checker']['objective'] = '0'
        self.highs[1]['claimed_status'] = 'time_limit'
        result = evaluate(self.manifest, self.lab[:3], self.highs[:3], self.campaign, panel=panel)
        self.assertTrue(result['panel_records_complete'])
        self.assertEqual(result['completed_correctly'], 2)
        self.assertEqual(result['strict_candidate_wins'], 1)
        self.assertFalse(result['completion_gate_passed'])
        self.assertFalse(result['speed_gate_passed'])
        with self.assertRaises(ValueError):
            evaluate(self.manifest, self.lab[:3], panel=[panel[0], panel[0]])

    def test_119_wins_and_ties_do_not_pass(self):
        for record in self.lab[119:]:
            record['wall_s'] = 2
        result = self.gate()
        self.assertTrue(result['completion_gate_passed'])
        self.assertFalse(result['passed'])
        self.assertEqual(result['strict_candidate_wins'], 119)
        self.assertEqual(result['ties'], 121)

    def test_suspension_invalidates_timing_for_either_solver(self):
        for runs in (self.lab, self.highs):
            runs[0]['sleep_affected'] = True
            result = self.gate()
            self.assertTrue(result['completion_gate_passed'])
            self.assertFalse(result['timing_comparable'])
            self.assertFalse(result['speed_gate_passed'])
            del runs[0]['sleep_affected']

    def test_missing_manifest_member_or_attempt_is_not_239_of_239(self):
        self.manifest['instances'].pop()
        self.lab.pop()
        self.assertFalse(self.gate()['completion_gate_passed'])
        self.manifest, self.lab, self.highs, self.campaign = fixture()
        self.lab.pop()
        result = self.gate()
        self.assertFalse(result['passed'])
        self.assertEqual(result['missing_candidate_instances'], ['m239'])

    def test_duplicates_do_not_hide_missing_models(self):
        self.lab[-1] = copy.deepcopy(self.lab[0])
        self.assertTrue(any('duplicate' in item for item in self.gate()['candidate_attempt_errors']))

    def test_mixed_candidate_binaries_cannot_complete_the_inventory(self):
        self.lab[17]['solver_sha256'] = 'e' * 64
        result = evaluate(self.manifest, self.lab)
        self.assertFalse(result['completion_gate_passed'])
        self.assertTrue(any('different solver binaries' in e for e in result['candidate_attempt_errors']))

    def test_unresolved_status_never_passes_even_with_reference_incumbent(self):
        for status in ['time_limit', 'node_limit', 'iteration_limit', 'unsupported', 'numerical_failure']:
            with self.subTest(status=status):
                self.lab[0]['claimed_status'] = status
                result = self.gate()
                self.assertFalse(result['completion_gate_passed'])
                self.assertEqual(result['completed_correctly'], 239)
                self.assertEqual(result['strict_candidate_wins'], 239)
                self.assertFalse(result['speed_gate_passed'])

    def test_invalid_unchecked_wrong_objective_and_weakened_checker_fail(self):
        changes = [('status', 'invalid'), ('checked_feasible', False), ('objective', '2'),
                   ('linear_tolerance', '.1'), ('integer_tolerance', '.1'),
                   ('objective_tolerance', '.1'), ('row_violations', 1),
                   ('max_row_violation', '1'), ('objective_error', '1')]
        for key, value in changes:
            with self.subTest(key=key):
                manifest, lab, highs, campaign = fixture()
                lab[0]['checker'][key] = value
                self.assertFalse(evaluate(manifest, lab, highs, campaign)['completion_gate_passed'])

    def test_reference_matching_point_requires_explicit_milp_checker_scope(self):
        checker = self.lab[0]['checker']
        for scope in (True, None, 0):
            if scope is None:
                del checker['relaxation']
            else:
                checker['relaxation'] = scope
            with self.subTest(scope=scope):
                result = self.gate()
                case = result['candidate_results']['m000']
                self.assertIn('original_MILP_checker_scope_required', case['reasons'])
                self.assertFalse(case['checked_feasible'])
                self.assertFalse(result['completion_gate_passed'])

    def test_unresolved_best_reference_is_not_upgraded_to_optimum(self):
        self.manifest['instances'][0]['reference_status'] = '=best='
        self.lab[0]['reference_status'] = self.highs[0]['reference_status'] = '=best='
        result = self.gate()
        self.assertFalse(result['passed'])
        self.assertEqual(result['unresolved_reference_instances'], ['m000'])

    def test_duplicate_model_hash_cannot_fill_manifest(self):
        self.manifest['instances'][1]['sha256'] = self.manifest['instances'][0]['sha256']
        self.assertFalse(self.gate()['completion_gate_passed'])

    def test_numerical_infeasibility_is_not_an_independent_proof(self):
        model = self.manifest['instances'][0]
        model.update(reference_status='=inf=', reference_objective='-')
        for records in [self.lab, self.highs]:
            records[0].update(reference_status='=inf=', reference_objective='-', claimed_status='infeasible',
                              checker_exit_code=None, checker={'status': 'unavailable', 'checked_feasible': False})
        self.assertFalse(self.gate()['completion_gate_passed'])
        model.update(validation_contract='reference_status_agreement', validation_basis='locked reference v37')
        result = self.gate()
        self.assertTrue(result['completion_gate_passed'])
        self.assertFalse(result['candidate_results']['m000']['independently_certified_infeasibility'])

    def test_header_only_infeasible_solution_is_an_absent_primal(self):
        model = self.manifest['instances'][0]
        model.update(reference_status='=inf=', reference_objective='-',
                     validation_contract='reference_status_agreement', validation_basis='locked v37')
        checker = dict(status='no_solution', checked_feasible=False, objective=None,
                       certificate_status='absent', certifies_infeasibility=False,
                       certifies_optimality=False)
        self.lab[0].update(reference_status='=inf=', reference_objective='-', claimed_status='infeasible',
                           checker_exit_code=4, checker=checker)
        self.assertTrue(self.gate()['completion_gate_passed'])
        for key, value in [('status', 'invalid'), ('status', 'error'), ('status', 'unsupported'),
                           ('checked_feasible', True), ('objective', '1'), ('certificate_status', 'invalid')]:
            with self.subTest(key=key, value=value):
                self.lab[0]['checker'] = dict(checker, **{key:value})
                self.assertFalse(self.gate()['completion_gate_passed'])
        self.lab[0]['checker'] = checker
        self.lab[0]['checker_exit_code'] = 2
        self.assertFalse(self.gate()['completion_gate_passed'])
        self.lab[0]['checker_exit_code'] = 4
        model['validation_contract'] = 'independent_original_model_infeasibility'
        self.assertFalse(self.gate()['completion_gate_passed'])

    def test_original_model_infeasibility_proof_can_satisfy_declared_contract(self):
        self.manifest['instances'][0].update(reference_status='=inf=', reference_objective='-')
        for records in [self.lab, self.highs]:
            records[0].update(reference_status='=inf=', reference_objective='-', claimed_status='infeasible',
                              checker={'status': 'valid', 'checked_feasible': False, 'relaxation': False, 'certificate_status': 'valid',
                                       'certifies_infeasibility': True})
        self.assertTrue(self.gate()['passed'])

    def test_no_historical_or_concurrent_engineering_speed_credit(self):
        for key, value in [('fair_timing_campaign', False), ('quiescent', False),
                           ('serial_solver_processes', False), ('coverage_exploration', True)]:
            with self.subTest(key=key):
                manifest, lab, highs, campaign = fixture()
                campaign[key] = value
                result = evaluate(manifest, lab, highs, campaign)
                self.assertTrue(result['completion_gate_passed'])
                self.assertEqual(result['strict_candidate_wins'], 0)
                self.assertFalse(result['passed'])
        self.assertFalse(evaluate(self.manifest, self.lab)['speed_gate_passed'])

    def test_mismatched_input_budget_campaign_hardware_configuration_rejected(self):
        for key, value in [('model_sha256', 'f'*64), ('campaign_id', 'historical'),
                           ('hardware_id', 'other-host'), ('time_limit_s', 20),
                           ('configuration_sha256', 'f'*64), ('solver_sha256', 'f'*64)]:
            with self.subTest(key=key):
                manifest, lab, highs, campaign = fixture()
                highs[0][key] = value
                result = evaluate(manifest, lab, highs, campaign)
                self.assertFalse(result['passed'])
                self.assertFalse(result['speed_comparisons'][0]['candidate_strictly_faster'])

    def test_timeout_baseline_never_becomes_candidate_speed_win(self):
        self.highs[0]['claimed_status'] = 'time_limit'
        result = self.gate()
        self.assertEqual(result['strict_candidate_wins'], 239)
        self.assertFalse(result['speed_comparisons'][0]['timing_eligible'])

    def test_relaxation_and_corrupt_json_rejected(self):
        self.lab[0]['relaxation'] = True
        self.assertFalse(self.gate()['completion_gate_passed'])
        for text in ['{"a":1,"a":2}', '{"x":NaN}', '{"x":Infinity}', '{"x":1e999}']:
            with self.assertRaises(ValueError):
                strict_json(text)


if __name__ == '__main__':
    unittest.main()
