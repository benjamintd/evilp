#!/usr/bin/env python3
"""Original-model GMI derivation contracts; no optimizer/engine imports."""
import argparse
import hashlib
from fractions import Fraction as Q
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

parser = argparse.ArgumentParser()
parser.add_argument('--checker', required=True, type=Path)
parser.add_argument('--primal-checker', type=Path)
args, rest = parser.parse_known_args()
if args.primal_checker is None:
    args.primal_checker = args.checker
counts = dict(certificate_calls=0, valid=0, invalid=0, unavailable=0,
              exact_mip_points=0, original_primal_checks=0)


def decimal(value):
    q = Q(value)
    sign, numerator, denominator = '-' if q < 0 else '', abs(q.numerator), q.denominator
    places = 0
    while denominator != 1:
        common = math.gcd(denominator, 10)
        if common == 1:
            raise ValueError(q)
        denominator //= common
        numerator *= 10 // common
        places += 1
    if not places:
        return sign + str(numerator)
    digits = str(numerator).zfill(places + 1)
    return sign + digits[:-places] + '.' + digits[-places:]


def mps(columns, rows, bounds, integers=(), duplicates=None):
    # rows: name, type, exact coefficient map, rhs, range|None.
    result = ['NAME INDEPENDENT', 'ROWS', ' N OBJ']
    result += [f' {kind} {name}' for name, kind, _, _, _ in rows]
    result += ['COLUMNS']
    for k, col in enumerate(columns):
        if col in integers:
            result += [f" M{k}A 'MARKER' 'INTORG'"]
        result += [f' {col} OBJ 0']
        for name, _, a, _, _ in rows:
            if (col, name) in (duplicates or {}):
                result += [f' {col} {name} {value}' for value in duplicates[col, name]]
            elif a.get(col, 0):
                result += [f' {col} {name} {decimal(a[col])}']
        if col in integers:
            result += [f" M{k}B 'MARKER' 'INTEND'"]
    result += ['RHS'] + [f' R {name} {decimal(rhs)}' for name, _, _, rhs, _ in rows]
    if any(r[-1] is not None for r in rows):
        result += ['RANGES'] + [f' R {name} {decimal(width)}'
                               for name, _, _, _, width in rows if width is not None]
    result += ['BOUNDS']
    for col in columns:
        lo, hi = bounds[col]
        result += [f' FR B {col}']
        if lo is not None:
            result += [f' LO B {col} {decimal(lo)}']
        if hi is not None:
            result += [f' UP B {col} {decimal(hi)}']
    return '\n'.join(result + ['ENDATA', ''])


def certificate(columns, box, beta, terms=(), equalities=(), fixed=(), published=None,
                rhs=Q(1, 2), basic='X'):
    lines = ['MILPLAB_GMI_1', 'SCOPE ORIGINAL_MIP_IN_BOX', f'BASIC {basic}', f'BETA {beta}']
    lines += [f"BOX {c} {'*' if box[c][0] is None else box[c][0]} "
              f"{'*' if box[c][1] is None else box[c][1]}" for c in columns]
    lines += [f'TERM {k} {a} {flag} {kind} {name} {side} {anchor}'
              for k, (a, flag, kind, name, side, anchor) in enumerate(terms)]
    lines += [f'EQ {name} {q}' for name, q in equalities]
    lines += [f'FIX {name} {q}' for name, q in fixed]
    lines += [f'PUB_RHS {rhs}']
    lines += [f'PUB {c} {(published or {}).get(c, 0)}' for c in columns]
    return '\n'.join(lines + ['END', ''])


def binary_cut(a, rhs, box):
    """Independent exact support compensation for nearest binary coefficients."""
    published = {c: Q(float(q)) for c, q in a.items()}
    correction = Q(0)
    for c, q in a.items():
        delta = published[c] - q
        correction += delta * (box[c][0] if delta > 0 else box[c][1] if delta < 0 else 0)
    maximum = rhs + correction
    rounded = float(maximum)
    if Q(rounded) > maximum:
        rounded = math.nextafter(rounded, -math.inf)
    return published, Q(rounded), correction


class GmiDerivation(unittest.TestCase):
    def check(self, model, cert, expected='valid', message=None):
        with tempfile.TemporaryDirectory(prefix='gmi-certificate-') as td:
            model_path, cert_path = Path(td)/'model.mps', Path(td)/'cut.txt'
            model_path.write_text(model)
            cert_path.write_bytes(cert.encode() if isinstance(cert, str) else cert)
            model_hash = hashlib.sha256(model_path.read_bytes()).hexdigest()
            certificate_hash = hashlib.sha256(cert_path.read_bytes()).hexdigest()
            process = subprocess.run([str(args.checker.resolve()), '--model', str(model_path),
                                      '--gmi', str(cert_path)], capture_output=True,
                                     text=True, timeout=25)
        self.assertEqual(process.stderr, '')
        data = json.loads(process.stdout)
        self.assertEqual(data['status'], expected, data)
        self.assertEqual(data['model_sha256'], model_hash)
        self.assertEqual(data['certificate_sha256'], certificate_hash)
        self.assertEqual(process.returncode, {'valid':0, 'invalid':1, 'unavailable':4,
                                             'unsupported':3}[expected], data)
        self.assertEqual(data['certifies_cut_validity'], expected == 'valid')
        for key in ('certifies_original_lp_infeasibility', 'certifies_original_mip_infeasibility',
                    'certifies_full_tree', 'recorded_box_proven_implied_by_original_model',
                    'checked_primal_witness', 'certifies_optimality'):
            self.assertFalse(data[key], (key, data))
        self.assertEqual(data['certificate_scope'], 'original_mip_in_recorded_box')
        if expected == 'valid':
            self.assertTrue(data['native_integrality_premises'])
            self.assertTrue(data['exact_tableau_identity'])
            self.assertGreaterEqual(Q(data['rounding_slack']), 0)
        if message:
            self.assertIn(message, data['message'])
        counts['certificate_calls'] += 1
        if expected in counts:
            counts[expected] += 1
        return data

    def primal(self, model, point, relax=False, expected=0):
        with tempfile.TemporaryDirectory(prefix='gmi-primal-control-') as td:
            model_path, sol = Path(td)/'model.mps', Path(td)/'point.sol'
            model_path.write_text(model)
            sol.write_text(''.join(f'{c} {q}\n' for c,q in point.items()))
            command = [str(args.primal_checker.resolve()), '--model', str(model_path),
                       '--solution', str(sol), '--format', 'rational', '--exact']
            if relax:
                command.append('--relax')
            process = subprocess.run(command, text=True, capture_output=True, timeout=20)
        self.assertEqual(process.stderr, '')
        self.assertEqual(process.returncode, expected, process.stdout)
        counts['original_primal_checks'] += 1
        return json.loads(process.stdout)

    def test_mode_arguments_and_strict_certificate_tokens(self):
        model = mps(['X'], [('E','E',{'X':Q(2)},Q(1),None)], {'X':(Q(0),Q(1))}, ['X'])
        cert = certificate(['X'], {'X':(Q(0),Q(1))}, Q(1,2), equalities=[('E',Q(1,2))])
        self.check(model, '# accepted comment\n\n' + cert)
        self.check(model, cert.replace('END\n', 'END $ ignored\n'), 'invalid')
        self.check(model, cert.replace('BETA 1/2', 'BETA 1/2 $ ignored'), 'invalid')
        self.check(model, cert.replace('PUB X 0', 'PUB X 0 # ignored'), 'invalid')
        with tempfile.TemporaryDirectory(prefix='gmi-mode-') as td:
            m, c = Path(td)/'model.mps', Path(td)/'cut.gmi'
            m.write_text(model)
            c.write_text(cert)
            base = [str(args.checker.resolve()), '--model', str(m), '--gmi', str(c)]
            variants = [ ['--inspect'], ['--solution', str(c)], ['--dual', str(c)],
                         ['--farkas', str(c)], ['--relax'], ['--exact'],
                         ['--format', 'simple'], ['--lin-tol','0'], ['--int-tol','0'],
                         ['--obj-tol','0'], ['--gmi',str(c)], ['--model',str(m)] ]
            for extra in variants:
                process = subprocess.run(base+extra, capture_output=True, text=True, timeout=20)
                self.assertEqual(process.returncode, 64, (extra,process.stdout,process.stderr))
                self.assertEqual(process.stdout, '')
                self.assertTrue(process.stderr)
            for command in ([str(args.checker.resolve()), '--gmi', str(c)],
                            [str(args.checker.resolve()), '--model', str(m), '--gmi']):
                process = subprocess.run(command, capture_output=True, text=True, timeout=20)
                self.assertEqual(process.returncode, 64)
            process = subprocess.run([str(args.checker.resolve()), '--help'], capture_output=True,
                                     text=True, timeout=20)
            self.assertEqual(process.returncode, 0)
            self.assertIn('--gmi', process.stdout)
            self.assertEqual(process.stderr, '')

    def test_lp_feasible_mip_in_box_infeasible_typing(self):
        model = mps(['X'], [('E','E',{'X':Q(2)},Q(1),None)], {'X':(Q(0),Q(1))}, ['X'])
        cert = certificate(['X'], {'X':(Q(0),Q(1))}, Q(1,2), equalities=[('E',Q(1,2))])
        data = self.check(model, cert)
        self.assertEqual(Q(data['exact_cut_rhs']), Q(1,2))
        self.assertEqual(data['exact_cut'], {'X':'0'})
        self.assertTrue(data['certifies_mip_in_box_infeasibility'])
        self.assertTrue(self.primal(model, {'X':Q(1,2)}, relax=True)['checked_feasible'])
        self.assertFalse(self.primal(model, {'X':Q(1,2)}, expected=1)['checked_feasible'])

    def test_mixed_integer_distances_negative_upper_and_shift(self):
        variants = [
            (Q(3,4), Q(1,2), Q(0), 'L', False, Q(3,4), Q(1,2)),
            (Q(3,4), Q(1,2), Q(0), 'L', True, Q(1,4), Q(1,2)),
            (Q(-3,4), Q(-1,2), Q(2), 'U', True, Q(-1,4), Q(0)),
            (Q(3,4), Q(-5,2), Q(4), 'L', True, Q(1,4), Q(3,2)),
        ]
        for a, beta, anchor, side, integer, gamma_original, expected_rhs in variants:
            sign = 1 if side == 'L' else -1
            source_a = sign*a
            source_rhs = beta + sign*a*anchor
            bounds = {'X':(None,None), 'Y':(anchor,anchor+8) if side == 'L' else (anchor-8,anchor)}
            box = {'X':(Q(-20),Q(20)), 'Y':bounds['Y']}
            model = mps(['X','Y'], [('E','E',{'X':Q(1),'Y':source_a},source_rhs,None)],
                        bounds, ['X','Y'] if integer else ['X'])
            cert = certificate(['X','Y'], box, beta,
                               [(a,'I' if integer else 'C','VAR','Y',side,anchor)], [('E',Q(1))],
                               published={'Y':gamma_original},rhs=expected_rhs)
            data = self.check(model, cert)
            self.assertEqual(Q(data['exact_cut']['Y']),gamma_original)
            self.assertEqual(Q(data['exact_cut_rhs']), expected_rhs)

    def test_exact_source_row_scaling_zero_ranges_and_other_valid_weakenings(self):
        box={'X':(Q(-10),Q(10)),'Y':(Q(0),Q(8))}
        for scale in (Q(1,1000000),Q(1),Q(32),Q(10**20)):
            model=mps(['X','Y'],[('E','E',{'X':scale,'Y':scale*Q(3,4)},scale/2,None)],
                      box,['X','Y'])
            cert=certificate(['X','Y'],box,Q(1,2),[(Q(3,4),'I','VAR','Y','L',Q(0))],
                             [('E',1/scale)],published={'Y':Q(1,4)},rhs=Q(1,2))
            self.check(model,cert)
            self.check(model,cert.replace('TERM 0 3/4','TERM 0 3/2'),'invalid','identity')
        model=mps(['X','Y'],[('E','G',{'X':Q(1),'Y':Q(3,4)},Q(1,2),Q(0))],
                  box,['X','Y'])
        cert=certificate(['X','Y'],box,Q(1,2),[(Q(3,4),'I','VAR','Y','L',Q(0))],
                         [('E',Q(1))],published={'Y':Q(1,4)},rhs=Q(1,2))
        self.check(model,cert)
        self.check(model,cert.replace('PUB Y 1/4','PUB Y 1/2'))
        self.check(model,cert.replace('PUB Y 1/4','PUB Y 1/8'),'invalid','published rhs')
        self.check(model,cert.replace('PUB Y 1/4','PUB Y 1/8').replace('PUB_RHS 1/2','PUB_RHS -1/2'))
        self.check(model,cert.replace('BOX Y 0 8','BOX Y 2 8'))  # Explicitly narrower scope.
        self.check(model.replace(' R E 0\n',' R E 1\n'),cert,'invalid','original equality')

    def test_row_slacks_ranges_and_bound_slacks(self):
        for kind, width, selected, anchor, a, beta, expected, rhs in [
                ('L',None,'U',Q(3),Q(1,2),Q(3,2),{'X':Q(-1)},Q(-1)),
                ('G',None,'L',Q(3),Q(-1,2),Q(3,2),{'X':Q(1),'Y':Q(1)},Q(2)),
                ('L',Q(2),'L',Q(1),Q(-1,2),Q(1,2),{'X':Q(1),'Y':Q(1)},Q(1))]:
            box = {'X':(Q(0),Q(10)), 'Y':(Q(0),Q(10))}
            model = mps(['X','Y'], [('R',kind,{'X':Q(2),'Y':Q(1)},Q(3),width)],box,['X'])
            cert = certificate(['X','Y'],box,beta,
                               [(a,'C','ROW','R',selected,anchor),
                                (Q(1,2),'C','VAR','Y','L',Q(0))],published=expected,rhs=rhs)
            self.check(model,cert)
            self.check(model,cert.replace(' C ROW ', ' I ROW '),'invalid','continuous')
        box = {'X':(Q(-10),Q(10)), 'Y':(Q(0),Q(2))}
        model = mps(['X','Y'], [('E','E',{'X':Q(1),'Y':Q(-3,4)},Q(-1),None)],box,['X'])
        cert = certificate(['X','Y'],box,Q(1,2),[(Q(3,4),'C','BOUND','Y','U',Q(2))],
                           [('E',Q(1))],published={'Y':Q(-3,4)},rhs=Q(-1))
        self.check(model,cert)
        self.check(model,cert.replace(' C BOUND ', ' I BOUND '),'invalid','native integer')

    def test_identity_modulo_fixed_box_is_explicitly_local(self):
        columns=['X','Y','F']
        bounds={'X':(None,None),'Y':(Q(0),Q(8)),'F':(Q(0),Q(4))}
        box={**bounds,'X':(Q(-20),Q(20)),'F':(Q(2),Q(2))}
        model=mps(columns,[('E','E',{'X':Q(1),'Y':Q(3,4),'F':Q(2)},Q(9,2),None)],bounds,['X'])
        cert=certificate(columns,box,Q(1,2),[(Q(3,4),'C','VAR','Y','L',Q(0))],
                         [('E',Q(1))],[('F',Q(-2))],{'Y':Q(3,4)},Q(1,2))
        self.check(model,cert)
        self.assertTrue(self.primal(model,{'X':Q(0),'Y':Q(0),'F':Q(9,4)})['checked_feasible'])
        self.check(model,cert.replace('FIX F -2','FIX F 2'),'invalid','identity')
        self.check(model,cert.replace('BOX F 2 2','BOX F 2 3'),'invalid','FIX')
        self.check(model,cert.replace('BOX F 2 2','BOX F 3 3'),'invalid','identity')

    def test_original_decimal_duplicates_and_safe_rounding(self):
        a, anchor, beta = Q(1,10), Q(-10**20), Q(1,2)
        bounds={'X':(None,None),'Y':(anchor,anchor+100)}
        box={**bounds,'X':(Q(-20),Q(20))}
        source_rhs=beta+a*anchor
        model=mps(['X','Y'],[('E','E',{'X':Q(1),'Y':a},source_rhs,None)],bounds,['X'],
                  {('Y','E'):['0.03','0.07']})
        exact_rhs=beta+a*anchor
        published,rhs,correction=binary_cut({'X':Q(0),'Y':a},exact_rhs,box)
        cert=certificate(['X','Y'],box,beta,[(a,'C','VAR','Y','L',anchor)],[('E',Q(1))],
                         published=published,rhs=rhs)
        data=self.check(model,cert)
        self.assertEqual(Q(data['rounding_error_lower']),correction)
        self.assertLess(correction,-500)
        self.check(model,cert.replace(f'PUB_RHS {rhs}',f'PUB_RHS {Q(float(exact_rhs))}'),
                   'invalid','published rhs')
        self.check(model,cert.replace(f'PUB Y {published["Y"]}','PUB Y 0.10000000000000001'),
                   'invalid','binary64')
        # Exact finite decimal serialization of the dyadic is accepted.
        self.check(model,cert.replace(f'PUB Y {published["Y"]}',f'PUB Y {decimal(published["Y"])}'))

    def test_tampered_integrality_native_anchors_coefficients_and_multipliers(self):
        box={'X':(Q(-10),Q(10)),'Y':(Q(0),Q(8))}
        model=mps(['X','Y'],[('E','E',{'X':Q(1),'Y':Q(3,4)},Q(1,2),None)],box,['X','Y'])
        cert=certificate(['X','Y'],box,Q(1,2),[(Q(3,4),'I','VAR','Y','L',Q(0))],
                         [('E',Q(1))],published={'Y':Q(1,4)},rhs=Q(1,2))
        self.check(model,cert)
        changes=[('BASIC X','BASIC NO_SUCH_COLUMN'),('BETA 1/2','BETA 1/4'),
                 ('TERM 0 3/4','TERM 0 1/2'),('EQ E 1','EQ E 0'),
                 ('VAR Y L 0','VAR Y L 1'),('VAR Y L 0','VAR Y U 0'),
                 ('VAR Y L 0','SCALED Y L 0'),('TERM 0','TERM 1'),
                 ('PUB Y 1/4','PUB Y 0'),('PUB Y 1/4','PUB Y -1/4'),
                 ('PUB X 0','PUB X 1'),('PUB_RHS 1/2',f'PUB_RHS {Q(math.nextafter(.5,math.inf))}'),
                 ('BOX Y 0 8','BOX Y 9 8')]
        for old,new in changes:
            self.check(model,cert.replace(old,new),'invalid')
        continuous=model.replace(" M1A 'MARKER' 'INTORG'\n",'').replace(" M1B 'MARKER' 'INTEND'\n",'')
        self.check(continuous,cert,'invalid','native integer')
        basic_continuous=model.replace(" M0A 'MARKER' 'INTORG'\n",'').replace(" M0B 'MARKER' 'INTEND'\n",'')
        self.check(basic_continuous,cert,'invalid','basic native')
        fractional_box={'X':box['X'],'Y':(Q(1,2),Q(8))}
        fractional_model=mps(['X','Y'],[('E','E',{'X':Q(1),'Y':Q(3,4)},Q(1,2),None)],
                             fractional_box,['X','Y'])
        fractional_cert=certificate(['X','Y'],fractional_box,Q(1,8),
                                   [(Q(3,4),'I','VAR','Y','L',Q(1,2))],[('E',Q(1))],
                                   published={'Y':Q(1,4)},rhs=Q(1,2))
        self.check(fractional_model,fractional_cert,'invalid','integer anchor')
        inequality=model.replace(' E E\n',' L E\n')
        self.check(inequality,cert,'invalid','original equality')

    def test_malformed_unknown_duplicates_and_missing_columns(self):
        model=mps(['X'],[('E','E',{'X':Q(2)},Q(1),None)],{'X':(Q(0),Q(1))},['X'])
        cert=certificate(['X'],{'X':(Q(0),Q(1))},Q(1,2),equalities=[('E',Q(1,2))])
        variants=[cert.replace('MILPLAB_GMI_1','MILPLAB_GMI_0'),
                  cert.replace('SCOPE ORIGINAL_MIP_IN_BOX','SCOPE ORIGINAL_LP'),
                  cert.replace('BOX X 0 1\n',''),cert.replace('PUB X 0\n',''),
                  cert.replace('BETA 1/2','BETA 1/0'),cert.replace('BETA 1/2','BETA nan'),
                  cert.replace('END\n',''),cert+'BETA 1/2\n',cert.replace('BOX X','BOX Z'),
                  cert.replace('PUB X','PUB Z'),cert.replace('EQ E','EQ OBJ'),
                  cert.replace('END','EQ E 1/2\nEND'),cert.replace('END','BOX X 0 1\nEND'),
                  cert.replace('END','PUB X 0\nEND'),cert.replace('END','BASIC X\nEND'),
                  cert.replace('END','FIX X 1\nEND'),cert.replace('END','WHAT 0\nEND'),
                  cert.encode().replace(b'BETA',b'BETA\x00')]
        for variant in variants:
            self.check(model,variant,'invalid')

    def test_integral_beta_and_unbounded_rounding_abstain(self):
        box={'X':(None,None),'Y':(Q(0),Q(8))}
        model=mps(['X','Y'],[('E','E',{'X':Q(1),'Y':Q(1,10)},Q(1),None)],box,['X'])
        cert=certificate(['X','Y'],box,Q(1),[(Q(1,10),'C','VAR','Y','L',Q(0))],
                         [('E',Q(1))],published={'Y':Q(.1)},rhs=Q(1,2))
        self.check(model,cert,'unavailable','integral native beta')
        model=model.replace(' R E 1\n',' R E 0.5\n')
        cert=cert.replace('BETA 1\n','BETA 1/2\n').replace('BOX Y 0 8','BOX Y * *')
        self.check(model,cert,'unavailable','missing recorded-box side')
        self.check(model,cert.replace(f'PUB Y {Q(.1)}','PUB Y 1/10'),'invalid','binary64')

    def test_resource_abstention_cleans_partial_certificate(self):
        model=mps(['X'],[('E','E',{'X':Q(2)},Q(1),None)],{'X':(Q(0),Q(1))},['X'])
        terms=[(Q(0),'C','VAR','X','L',Q(0))]*4097
        cert=certificate(['X'],{'X':(Q(0),Q(1))},Q(1,2),terms,[('E',Q(1,2))])
        self.check(model,cert,'unavailable','term capacity')

    def test_generated_exact_native_tableaux_and_fraction_points(self):
        rng=random.Random(77309)
        for k in range(72):
            a=Q(rng.randint(-15,15),rng.choice([2,4,5,8,10]))
            beta=Q(rng.randint(-15,15),rng.choice([4,5,8]))
            if beta.denominator==1:
                beta+=Q(1,4)
            anchor=Q(rng.randint(-4,4))
            integer=bool(k%2)
            f=beta-(beta.numerator//beta.denominator)
            if integer:
                fa=a-(a.numerator//a.denominator)
                gamma=fa if fa<=f else f*(1-fa)/(1-f)
            else:
                gamma=a if a>=0 else -f*a/(1-f)
            box={'X':(Q(-100),Q(100)),'Y':(anchor,anchor+8)}
            source_rhs=beta+a*anchor
            model=mps(['X','Y'],[('E','E',{'X':Q(1),'Y':a},source_rhs,None)],
                      {'X':(None,None),'Y':box['Y']},['X','Y'] if integer else ['X'])
            exact={'X':Q(0),'Y':gamma}
            exact_rhs=f+gamma*anchor
            published,rhs,correction=binary_cut(exact,exact_rhs,box)
            cert=certificate(['X','Y'],box,beta,[(a,'I' if integer else 'C','VAR','Y','L',anchor)],
                             [('E',Q(1))],published=published,rhs=rhs)
            data=self.check(model,cert)
            self.assertEqual({c:Q(q) for c,q in data['exact_cut'].items()},exact)
            self.assertEqual(Q(data['exact_cut_rhs']),exact_rhs)
            self.assertEqual(Q(data['rounding_error_lower']),correction)
            for x in map(Q,range(-40,41)):
                if not a:
                    continue
                y=(source_rhs-x)/a
                if not anchor<=y<=anchor+8 or (integer and y.denominator!=1):
                    continue
                self.assertEqual(x+a*y,source_rhs)
                self.assertGreaterEqual(gamma*y,exact_rhs)
                self.assertGreaterEqual(published['Y']*y,rhs)
                counts['exact_mip_points']+=1
        self.assertGreater(counts['exact_mip_points'],100)


if __name__=='__main__':
    result=unittest.main(argv=[sys.argv[0]]+rest,exit=False)
    print(json.dumps({'validation_counts':counts},sort_keys=True))
    sys.exit(not result.result.wasSuccessful())
