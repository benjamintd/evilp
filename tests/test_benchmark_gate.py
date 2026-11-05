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
                'fair_timing_campaign': True, 'serial_solver_processes': True, 'quiescent': True,
                'coverage_exploration': False, 'measurement': MEASUREMENT, 'budget': budget,
                'engines': {name: {'solver_sha256': char*64, 'configuration_sha256': config*64}
                            for name, char, config in [('lab', 'a', 'b'), ('highs', 'c', 'd')]}}
    def runs(solver):
        return [{'instance': model['id'], 'solver': solver, 'model_sha256': model['sha256'],
                 'mode': 'mip', 'relaxation': False, 'claimed_status': 'optimal', 'exit_code': 0,
                 'signal': 0, 'watchdog_timeout': False, 'status_conflict': False,
                 'reference_status': '=opt=', 'reference_objective': '1', 'checker_exit_code': 0,
                 'checker': {'status': 'valid', 'checked_feasible': True, 'objective': '1',
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

    def test_all_240_complete_and_120_strict_wins_pass(self):
        for record in self.lab[120:]:
            record['wall_s'] = 3
        result = self.gate()
        self.assertTrue(result['passed'])
        self.assertEqual(result['strict_candidate_wins'], 120)
        self.assertEqual(result['speed_denominator'], 240)
        self.assertFalse(result['candidate_results']['m000']['independently_certified_optimality'])

    def test_119_wins_and_ties_do_not_pass(self):
        for record in self.lab[119:]:
            record['wall_s'] = 2
        result = self.gate()
        self.assertTrue(result['completion_gate_passed'])
        self.assertFalse(result['passed'])
        self.assertEqual(result['strict_candidate_wins'], 119)
        self.assertEqual(result['ties'], 121)

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

    def test_original_model_infeasibility_proof_can_satisfy_declared_contract(self):
        self.manifest['instances'][0].update(reference_status='=inf=', reference_objective='-')
        for records in [self.lab, self.highs]:
            records[0].update(reference_status='=inf=', reference_objective='-', claimed_status='infeasible',
                              checker={'status': 'valid', 'checked_feasible': False, 'certificate_status': 'valid',
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
