#!/usr/bin/env python3
"""Independent arithmetic oracle and adversarial parser/solution contracts."""
import fractions, gzip, json, pathlib, random, subprocess, sys, tempfile, unittest
BIN=str(pathlib.Path(sys.argv.pop(1)).resolve())
FIX=pathlib.Path(__file__).parent/'fixtures'
class Checker(unittest.TestCase):
    def call(self,model,solution=None,fmt='simple',extra=()):
        args=[BIN,'--model',str(model)]
        args+=['--inspect'] if solution is None else ['--solution',str(solution),'--format',fmt]
        p=subprocess.run(args+list(extra),capture_output=True,text=True)
        self.assertIn(p.returncode,range(5),p.stderr)
        return p.returncode,json.loads(p.stdout)
    def check_text(self,sol,expected=0,model=None,fmt='simple',extra=()):
        with tempfile.TemporaryDirectory() as td:
            p=pathlib.Path(td)/'x.sol';p.write_text(sol)
            code,result=self.call(model or FIX/'tiny.mps',p,fmt,extra)
            self.assertEqual(code,expected,result)
            return result
    def test_exact_decimal_offset_and_sparse_zero(self):
        r=self.check_text((FIX/'valid.sol').read_text(),extra=['--exact'])
        self.assertEqual(r['objective'],'803/100')
        self.assertTrue(r['checked_feasible'])
        self.assertFalse(r['certifies_optimality'])
    def test_row_bound_integer_objective_rejections(self):
        for sol,key in [('x 1\ny 1\nz .3\n','row_violations'),('x 2\nz .3\n','bound_violations'),('x .5\ny 1\nz .3\n','integrality_violations')]:
            self.assertGreater(self.check_text(sol,1)[key],0)
        self.check_text('=obj= 8.04\nx 1\nz .3\n',1)
    def test_malformed_nonfinite_duplicate_unknown(self):
        for sol in ['','x nan\n','x inf\n','x 1e100000\n','x 1\nx 1\n','unknown 1\n','x garbage\n']:
            self.check_text(sol,2)
    def test_no_solution_is_not_a_certificate(self):
        for fmt,sol in [('simple','=infeas=\n'),('scip','solution status: infeasible\n'),('scip','solution status: time limit reached\nno solution available\n'),('cbc','Infeasible - objective value 1\n'),('highs','Model status\nInfeasible\n')]:
            r=self.check_text(sol,4,fmt=fmt);self.assertFalse(r['checked_feasible'])
    def test_adapter_formats(self):
        examples={
          'scip':'solution status: optimal solution found\nobjective value: 8.03\nx 1 (obj:3)\nz .3 (obj:0.1)\n',
          'cbc':'Optimal - objective value 8.03000000\n 0 x 1 3\n 1 y 0 2\n 2 z .3 .1\n',
          'highs':'Model status\nOptimal\n# Primal solution values\nFeasible\nObjective 8.03\n# Columns 3\nx 1\ny 0\nz .3\n# Rows 2\nCAP 2\nFIXZ .3\n# Dual solution values\nNone\n'}
        for fmt,sol in examples.items():self.check_text(sol,fmt=fmt)
    def test_tolerance_is_exact_and_explicit(self):
        self.check_text('x 1\nz .30001\n',extra=['--lin-tol','0.00001'])
        self.check_text('x 1\nz .3000100000000000000001\n',1,extra=['--lin-tol','0.00001'])
        self.check_text('x 1\nz .3000000000000000000001\n',1,extra=['--exact'])
    def test_gzip_and_d_exponents(self):
        with tempfile.TemporaryDirectory() as td:
            m=pathlib.Path(td)/'model.mps.gz'
            with gzip.open(m,'wt') as f:f.write((FIX/'tiny.mps').read_text().replace('0.3','3D-1'))
            self.check_text('=obj= 8.03\nx 1\nz .3\n',model=m,extra=['--exact'])
    def test_ranges_all_senses(self):
        for sense,ran,valid,invalid in [('L','-2','4','2'),('G','-2','6','8'),('E','2','6','4'),('E','-2','4','6')]:
            model=f'NAME R\nROWS\n N O\n {sense} R\nCOLUMNS\n x O 0 R 1\nRHS\n b R 5\nRANGES\n g R {ran}\nBOUNDS\n FR b x\nENDATA\n'
            with tempfile.TemporaryDirectory() as td:
                m=pathlib.Path(td)/'r.mps';m.write_text(model)
                self.check_text('x '+valid+'\n',model=m,extra=['--exact'])
                self.check_text('x '+invalid+'\n',1,model=m,extra=['--exact'])
    def test_exact_cancellation_beyond_binary_float_precision(self):
        model='NAME CANCEL\nROWS\n N O\n E R\nCOLUMNS\n x O 0 R 1e30\n z R 1\n y R -1e30\nRHS\n b R 0\nBOUNDS\n FR b x\n FR b y\n FR b z\nENDATA\n'
        with tempfile.TemporaryDirectory() as td:
            m=pathlib.Path(td)/'cancel.mps';m.write_text(model)
            self.check_text('x 1\ny 1\nz 0\n',model=m,extra=['--exact'])
            r=self.check_text('x 1\ny 1\nz 0.00001\n',1,model=m,extra=['--exact'])
            self.assertEqual(r['max_row_violation'],'1/100000')
    def test_integer_defaults_and_explicit_bounds(self):
        for bound,value,code in [('', '2',1),(' LO b x 2\n','2',0),(' LI b x 2\n UI b x 4\n','3',0),(' BV b x\n','0.5',1),(' BV b x 1\n','1',0)]:
            model='NAME I\nROWS\n N O\nCOLUMNS\n M0 \'MARKER\' \'INTORG\'\n x O 1\n M1 \'MARKER\' \'INTEND\'\nBOUNDS\n'+bound+'ENDATA\n'
            with tempfile.TemporaryDirectory() as td:
                m=pathlib.Path(td)/'i.mps';m.write_text(model);self.check_text('x '+value+'\n',code,model=m,extra=['--exact'])
    def test_conflicting_marker_declarations_fail_closed(self):
        text="NAME BAD\nROWS\n N O\n E R\nCOLUMNS\n x O 1\n M0 'MARKER' 'INTORG'\n x R 1\n M1 'MARKER' 'INTEND'\nRHS\n b R 0.5\nENDATA\n"
        with tempfile.TemporaryDirectory() as td:
            m=pathlib.Path(td)/'conflict.mps';m.write_text(text)
            code,r=self.call(m);self.assertEqual(code,2,r)
    def test_fail_closed_model_extensions_and_truncation(self):
        original=(FIX/'tiny.mps').read_text()
        with tempfile.TemporaryDirectory() as td:
            for extension in [' SC bounds z 1\n',' SI bounds z 1\n']:
                m=pathlib.Path(td)/'bad.mps';m.write_text(original.replace(' FR bounds z\n',extension))
                self.assertEqual(self.call(m)[0],3)
            m.write_text(original.replace('ENDATA','QMATRIX\n x x 2\nENDATA'))
            self.assertEqual(self.call(m)[0],3)
            m.write_text(original.replace('ENDATA',''));self.assertEqual(self.call(m)[0],2)
    def test_randomized_exact_linear_oracle(self):
        rng=random.Random(20261105)
        F=fractions.Fraction
        with tempfile.TemporaryDirectory() as td:
            m=pathlib.Path(td)/'r.mps';s=pathlib.Path(td)/'r.sol'
            for _ in range(120):
                a=[rng.randint(-15,15) for _ in range(4)];x=[rng.randint(-20,20) for _ in range(4)]
                rhs=rng.randint(-100,100);sense=rng.choice('LEG')
                lhs=sum(F(ai,10)*xi for ai,xi in zip(a,x))
                valid={'L':lhs<=rhs,'E':lhs==rhs,'G':lhs>=rhs}[sense]
                coeff=lambda n:f'{n/10:.1f}'
                m.write_text('NAME RAND\nROWS\n N O\n '+sense+' R\nCOLUMNS\n'+''.join(f' x{i} O {coeff(a[i])} R {coeff(a[i])}\n' for i in range(4))+f'RHS\n b R {rhs}\nBOUNDS\n'+''.join(f' FR b x{i}\n' for i in range(4))+'ENDATA\n')
                s.write_text(''.join(f'x{i} {x[i]}\n' for i in range(4)))
                code,result=self.call(m,s,extra=['--exact'])
                self.assertEqual(code,0 if valid else 1,result)
                self.assertEqual(F(result['objective']),lhs)
unittest.main()
