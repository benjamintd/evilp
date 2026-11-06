#!/usr/bin/env python3
"""Independently replay engine-format root cuts and scoped MIP-box proofs."""
import argparse
import copy
from fractions import Fraction as Q
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

parser = argparse.ArgumentParser()
parser.add_argument('--checker', type=Path, required=True)
args, rest = parser.parse_known_args()
CHECKER = str(args.checker.resolve())
REPLAY = Path(__file__).resolve().parents[1] / 'scripts/check_trace_certificates.py'
COUNTS = dict(replays=0, accepted=0, rejected=0, exact_cut_checks=0)

# Retained output of the native exporter -> exact GMI -> publication pipeline:
# integer X, 2X<=3, 0<=X<=16 produces -X>=-1. This is not a stored answer
# used by the solver; the independent checker re-derives its complete proof.
SCALAR = ("NAME S\nROWS\n N OBJ\n L R\nCOLUMNS\n M0 'MARKER' 'INTORG'\n"
          " X OBJ -1 R 2\n M1 'MARKER' 'INTEND'\nRHS\n RHS R 3\n"
          "BOUNDS\n UP B X 16\nENDATA\n")
SCALAR_CERT = ('MILPLAB_GMI_1\nSCOPE ORIGINAL_MIP_IN_BOX\nBASIC X\nBETA 3/2\n'
               'BOX X 0 16\nTERM 0 1/2 C ROW R U 3\nPUB_RHS -1\nPUB X -1\nEND\n')


def root_cut(certificate=SCALAR_CERT, name='CUT0'):
    return dict(event='root_gmi_cut', cut_id=name, scope='ORIGINAL_MIP_IN_BOX',
                certificate=certificate)


def proof(name='CUT0', lower='5/4', upper='3/2', margin='1/4'):
    return dict(event='mip_in_box_infeasibility_certificate', lp_solve=2,
                cut_ids=[name], bounds=[dict(column='X', lower=lower, upper=upper)],
                multipliers={name:'1'}, contradiction_margin=margin)


class CutTrace(unittest.TestCase):
    def replay(self, directory, model, events, expected=0, compressed=False, extra=()):
        directory.mkdir()
        source = directory / ('original.mps.gz' if compressed else 'original.mps')
        source.write_bytes(gzip.compress(model.encode()) if compressed else model.encode())
        trace, output = directory/'trace.jsonl', directory/'replayed'
        trace.write_text(events if isinstance(events,str) else
                         ''.join(json.dumps(e)+'\n' for e in events))
        command = [sys.executable,str(REPLAY),'--model',str(source),'--trace',str(trace),
                   '--checker',CHECKER,'--out',str(output),*extra]
        process = subprocess.run(command,capture_output=True,text=True,timeout=50)
        self.assertEqual(process.returncode,expected,(process.stdout,process.stderr))
        self.assertEqual(process.stderr,'')
        summary = json.loads(process.stdout)
        self.assertEqual(summary,json.loads((output/'summary.json').read_text()))
        self.assertFalse(summary['certifies_tree_infeasibility'])
        self.assertFalse(summary['certifies_original_lp_infeasibility'])
        self.assertFalse(summary['certifies_original_mip_infeasibility'])
        self.assertEqual(summary['all_certificates_valid'],expected==0)
        if summary['cut_count'] or summary['mip_in_box_proof_count']:
            self.assertFalse(summary['all_lp_certificates_valid'])
        for record in summary['proofs']:
            self.assertFalse(record['certifies_full_tree'])
            self.assertFalse(record['certifies_original_lp_infeasibility'])
            self.assertFalse(record['certifies_original_mip_infeasibility'])
        COUNTS['replays']+=1
        COUNTS['accepted' if expected==0 else 'rejected']+=1
        COUNTS['exact_cut_checks']+=summary['cut_count']
        return summary,output

    def test_actual_engine_scalar_cut_and_cut_dependent_proof(self):
        with tempfile.TemporaryDirectory(prefix='cut-trace-') as td:
            summary,out=self.replay(Path(td)/'run',SCALAR,[root_cut(),proof()],compressed=True)
            self.assertEqual(summary['valid_cut_count'],1)
            self.assertEqual(summary['valid_proof_count'],1)
            self.assertEqual(summary['mip_in_box_proof_count'],1)
            self.assertEqual(summary['cuts'][0]['checker']['exact_cut'],{'X':'-1'})
            p=summary['proofs'][0]
            self.assertEqual(Q(p['checked_margin']),Q(1,4))
            self.assertEqual(p['certificate_scope'],'original_mip_in_recorded_box')
            model=(out/p['directory']/'model.mps').read_text()
            # The new X entry stays inside the existing integer marker pair.
            self.assertIn(" X OBJ -1 R 2\n X CUT0 -1\n M1 'MARKER' 'INTEND'",model)
            self.assertIn(' RHS CUT0 -1\n',model)

    def test_retained_actual_branch_search_trace_has_two_scoped_contradictions(self):
        # Generated seed2026110602, second six-binary/four-row instance. Actual
        # search emits three cuts and two Farkas certificates; the original
        # MIP also has the independently checked feasible point (1,1,0,0,0,0).
        model = 'NAME GENERATED\nOBJSENSE\n MAX\nROWS\n N O\n L R0\n L R1\n L R2\n L R3\nCOLUMNS\n X0 O 6\n X0 R0 2\n X0 R1 -1\n X0 R2 5\n X0 R3 -3\n X1 O 3\n X1 R0 -2\n X1 R1 -4\n X1 R2 3\n X1 R3 -2\n X2 O 3\n X2 R0 1\n X2 R1 4\n X2 R2 -1\n X2 R3 2\n X3 O 4\n X3 R0 1\n X3 R1 3\n X3 R2 1\n X3 R3 -2\n X4 O 2\n X4 R0 4\n X4 R1 -2\n X4 R2 2\n X4 R3 2\n X5 O 7\n X5 R0 3\n X5 R1 4\n X5 R3 4\nRHS\n rhs R0 4\n rhs R1 -3\n rhs R2 9\n rhs R3 -2\nBOUNDS\n BV box X0\n BV box X1\n BV box X2\n BV box X3\n BV box X4\n BV box X5\nENDATA\n'
        events = json.loads(r'''[
  {
    "event": "root_gmi_cut",
    "cut_id": "__milplab_gmi_0",
    "scope": "ORIGINAL_MIP_IN_BOX",
    "certificate": "MILPLAB_GMI_1\nSCOPE ORIGINAL_MIP_IN_BOX\nBASIC X3\nBETA 1/7\nBOX X0 0 1\nBOX X1 0 1\nBOX X2 0 1\nBOX X3 0 1\nBOX X4 0 1\nBOX X5 0 1\nTERM 0 1/7 C ROW R1 U -3\nTERM 1 2/7 C ROW R2 U 9\nTERM 2 -1/7 C ROW R3 U -2\nTERM 3 -12/7 C BOUND X0 U 1\nTERM 4 -4/7 C BOUND X1 U 1\nPUB_RHS -5254199565265579/2251799813685248\nPUB X0 -3/2\nPUB X1 -6004799503160661/18014398509481984\nPUB X2 -6004799503160661/18014398509481984\nPUB X3 -6004799503160661/9007199254740992\nPUB X4 -6004799503160661/18014398509481984\nPUB X5 -6004799503160661/9007199254740992\nEND\n"
  },
  {
    "event": "root_gmi_cut",
    "cut_id": "__milplab_gmi_1",
    "scope": "ORIGINAL_MIP_IN_BOX",
    "certificate": "MILPLAB_GMI_1\nSCOPE ORIGINAL_MIP_IN_BOX\nBASIC X4\nBETA 3/7\nBOX X0 0 1\nBOX X1 0 1\nBOX X2 0 1\nBOX X3 0 1\nBOX X4 0 1\nBOX X5 0 1\nTERM 0 -1/2 I VAR X2 L 0\nTERM 1 -1/14 C ROW R1 U -3\nTERM 2 5/14 C ROW R2 U 9\nTERM 3 1/14 C ROW R3 U -2\nTERM 4 -23/14 C BOUND X0 U 1\nTERM 5 -17/14 C BOUND X1 U 1\nPUB_RHS -37/8\nPUB X0 -11/4\nPUB X1 -13/8\nPUB X2 3/8\nPUB X3 -3/8\nPUB X4 -3/4\nPUB X5 -1/2\nEND\n"
  },
  {
    "event": "root_gmi_cut",
    "cut_id": "__milplab_gmi_2",
    "scope": "ORIGINAL_MIP_IN_BOX",
    "certificate": "MILPLAB_GMI_1\nSCOPE ORIGINAL_MIP_IN_BOX\nBASIC X5\nBETA 17/28\nBOX X0 0 1\nBOX X1 0 1\nBOX X2 0 1\nBOX X3 0 1\nBOX X4 0 1\nBOX X5 0 1\nTERM 0 3/4 I VAR X2 L 0\nTERM 1 3/28 C ROW R1 U -3\nTERM 2 -1/28 C ROW R2 U 9\nTERM 3 1/7 C ROW R3 U -2\nTERM 4 5/7 C BOUND X0 U 1\nTERM 5 23/28 C BOUND X1 U 1\nPUB_RHS -7369526662969903/9007199254740992\nPUB X0 -8188362958855447/18014398509481984\nPUB X1 -1228254443828317/4503599627370496\nPUB X2 -1228254443828317/4503599627370496\nPUB X3 -3275345183542179/36028797018963968\nPUB X4 -3275345183542179/18014398509481984\nPUB X5 -1\nEND\n"
  },
  {
    "event": "mip_in_box_infeasibility_certificate",
    "lp_solve": 9,
    "scope": "Original MIP with recorded bounds and GMI cuts; not original LP or a tree proof",
    "cut_ids": [
      "__milplab_gmi_0",
      "__milplab_gmi_1",
      "__milplab_gmi_2"
    ],
    "bounds": [
      {
        "column": "X0",
        "lower": "0",
        "upper": "1"
      },
      {
        "column": "X1",
        "lower": "0",
        "upper": "1"
      },
      {
        "column": "X2",
        "lower": "1",
        "upper": "1"
      },
      {
        "column": "X3",
        "lower": "0",
        "upper": "0"
      },
      {
        "column": "X4",
        "lower": "0",
        "upper": "0"
      },
      {
        "column": "X5",
        "lower": "0",
        "upper": "1"
      }
    ],
    "multipliers": {
      "R0": "0",
      "R1": "-1/4",
      "R2": "0",
      "R3": "-1/4",
      "__milplab_gmi_0": "0",
      "__milplab_gmi_1": "0",
      "__milplab_gmi_2": "18014398509481984/8188362958855447"
    },
    "contradiction_margin": "21289743693024155/32753451835421788"
  },
  {
    "event": "mip_in_box_infeasibility_certificate",
    "lp_solve": 10,
    "scope": "Original MIP with recorded bounds and GMI cuts; not original LP or a tree proof",
    "cut_ids": [
      "__milplab_gmi_0",
      "__milplab_gmi_1",
      "__milplab_gmi_2"
    ],
    "bounds": [
      {
        "column": "X0",
        "lower": "0",
        "upper": "1"
      },
      {
        "column": "X1",
        "lower": "0",
        "upper": "1"
      },
      {
        "column": "X2",
        "lower": "0",
        "upper": "1"
      },
      {
        "column": "X3",
        "lower": "1",
        "upper": "1"
      },
      {
        "column": "X4",
        "lower": "0",
        "upper": "0"
      },
      {
        "column": "X5",
        "lower": "0",
        "upper": "1"
      }
    ],
    "multipliers": {
      "R0": "0",
      "R1": "-1/4",
      "R2": "0",
      "R3": "0",
      "__milplab_gmi_0": "1/6",
      "__milplab_gmi_1": "0",
      "__milplab_gmi_2": "0"
    },
    "contradiction_margin": "30023997515803303/108086391056891904"
  }
]''')
        with tempfile.TemporaryDirectory() as td:
            s,out=self.replay(Path(td)/'run',model,events)
            self.assertEqual(s['valid_cut_count'],3)
            self.assertEqual(s['valid_proof_count'],2)
            self.assertEqual(s['mip_in_box_proof_count'],2)
            point=Path(td)/'point.sol'
            point.write_text('X0 1\nX1 1\nX2 0\nX3 0\nX4 0\nX5 0\n')
            p=subprocess.run([CHECKER,'--model',str(out.parent/'original.mps'),'--solution',str(point),'--exact'],
                             capture_output=True,text=True,timeout=20)
            self.assertEqual(p.returncode,0,p.stdout)
            self.assertTrue(json.loads(p.stdout)['checked_feasible'])

    def test_cut_only_trace_checks_derivation_without_infeasibility_claim(self):
        with tempfile.TemporaryDirectory() as td:
            s,_=self.replay(Path(td)/'run',SCALAR,[root_cut()])
            self.assertEqual(s['cut_count'],1)
            self.assertEqual(s['proof_count'],0)
            self.assertFalse(s['cuts'][0]['checker']['certifies_mip_in_box_infeasibility'])

    def test_parity_constant_cut_is_mip_only_while_original_lp_is_feasible(self):
        model=('NAME P\nROWS\n N O\n E E\nCOLUMNS\n M0 \'MARKER\' \'INTORG\'\n'
               ' X O 0 E 2\n M1 \'MARKER\' \'INTEND\'\nRHS\n b E 1\n'
               'BOUNDS\n UP b X 1\nENDATA\n')
        cert=('MILPLAB_GMI_1\nSCOPE ORIGINAL_MIP_IN_BOX\nBASIC X\nBETA 1/2\n'
              'BOX X 0 1\nEQ E 1/2\nPUB_RHS 1/2\nPUB X 0\nEND\n')
        event=proof(lower='0',upper='1',margin='1/2')
        with tempfile.TemporaryDirectory() as td:
            s,out=self.replay(Path(td)/'run',model,[root_cut(cert),event])
            self.assertTrue(s['cuts'][0]['checker']['certifies_mip_in_box_infeasibility'])
            sol=Path(td)/'lp.sol';sol.write_text('X 1/2\n')
            p=subprocess.run([CHECKER,'--model',str(out.parent/'original.mps'),'--solution',str(sol),
                              '--format','rational','--relax','--exact'],capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stdout)
            self.assertTrue(json.loads(p.stdout)['checked_feasible'])
            self.assertEqual(Q(s['proofs'][0]['checked_margin']),Q(1,2))

    def test_publication_preserves_full_binary_dyadics_and_existing_rhs_vector(self):
        coefficient=-Q(1)-Q(1,2**52)
        rhs=-Q(1)-Q(1,2**48)
        cert=SCALAR_CERT.replace('PUB X -1',f'PUB X {coefficient}').replace('PUB_RHS -1',f'PUB_RHS {rhs}')
        event=proof(margin=str(rhs-coefficient*Q(5,4)))
        with tempfile.TemporaryDirectory() as td:
            s,out=self.replay(Path(td)/'run',SCALAR,[root_cut(cert),event])
            model=(out/s['proofs'][0]['directory']/'model.mps').read_text()
            self.assertIn(' X CUT0 -1.0000000000000002220446049250313080847263336181640625\n',model)
            self.assertIn(' RHS CUT0 -1.000000000000003552713678800500929355621337890625\n',model)
            self.assertEqual(Q(s['proofs'][0]['checked_margin']),Q(event['contradiction_margin']))

    def test_source_decimals_duplicates_ranges_continuations_and_bound_only_columns(self):
        model=("NAME S\nOBJSENSE\n MAX\nROWS\n N OBJ\n L R\nCOLUMNS\n"
               " M0 'MARKER' 'INTORG'\n X OBJ -1 R .7\n R 1.3\n"
               " M1 'MARKER' 'INTEND'\nRHS\n RHS OBJ 7 R 3\nRANGES\n range R 1\n"
               "BOUNDS\n UP B X 16\n FR B Z\nENDATA\n")
        cert=SCALAR_CERT.replace('BOX X 0 16','BOX X 0 16\nBOX Z * *').replace('PUB X -1','PUB X -1\nPUB Z 0')
        event=proof();event['bounds'].append(dict(column='Z',lower=None,upper=None))
        with tempfile.TemporaryDirectory() as td:
            s,out=self.replay(Path(td)/'run',model,[root_cut(cert),event],compressed=True)
            rewritten=(out/s['proofs'][0]['directory']/'model.mps').read_text()
            self.assertIn(' X OBJ -1 R .7\n R 1.3\n X CUT0 -1\n',rewritten)
            self.assertIn('RHS\n RHS OBJ 7 R 3\n RHS CUT0 -1\nRANGES\n range R 1\n',rewritten)
            self.assertIn(' FR replay Z\n',rewritten)

    def test_original_without_rhs_gets_one_exact_rhs_section(self):
        model=("NAME NORHS\nROWS\n N O\n E E\nCOLUMNS\n M0 'MARKER' 'INTORG'\n"
               " X O 0 E 2\n M1 'MARKER' 'INTEND'\n Y O 0 E 1\n"
               "BOUNDS\n FR b X\n LO b Y -1.5\nENDATA\n")
        cert=('MILPLAB_GMI_1\nSCOPE ORIGINAL_MIP_IN_BOX\nBASIC X\nBETA 3/4\n'
              'BOX X * *\nBOX Y -3/2 *\nTERM 0 1/2 C VAR Y L -3/2\nEQ E 1/2\n'
              'PUB_RHS 0\nPUB X 0\nPUB Y 1/2\nEND\n')
        event=proof(margin='1/4')
        event['bounds']=[dict(column='X',lower='0',upper='1'),dict(column='Y',lower='-1',upper='-1/2')]
        with tempfile.TemporaryDirectory() as td:
            s,out=self.replay(Path(td)/'run',model,[root_cut(cert),event])
            rewritten=(out/s['proofs'][0]['directory']/'model.mps').read_text()
            self.assertEqual(rewritten.splitlines().count('RHS'),1)
            self.assertIn('RHS\n replaycuts CUT0 0\n',rewritten)

    def test_multiple_proved_cut_rows_are_inserted_and_named_independently(self):
        model=("NAME TWO\nROWS\n N OBJ\n L R\n L S\nCOLUMNS\n M0 'MARKER' 'INTORG'\n"
               " X OBJ -1 R 2\n Y OBJ -1 S 2\n M1 'MARKER' 'INTEND'\nRHS\n RHS R 3 S 3\n"
               "BOUNDS\n UP B X 16\n UP B Y 16\nENDATA\n")
        first=SCALAR_CERT.replace('BOX X 0 16','BOX X 0 16\nBOX Y 0 16').replace('PUB X -1','PUB X -1\nPUB Y 0')
        second=(first.replace('BASIC X','BASIC Y').replace('ROW R U','ROW S U')
                .replace('PUB X -1','PUB X 0').replace('PUB Y 0','PUB Y -1'))
        event=proof(margin='1/2')
        event['cut_ids']=['CUT0','CUT1']
        event['bounds'].append(dict(column='Y',lower='5/4',upper='3/2'))
        event['multipliers']['CUT1']='1'
        with tempfile.TemporaryDirectory() as td:
            good=[root_cut(first),root_cut(second,'CUT1'),event]
            s,out=self.replay(Path(td)/'good',model,good)
            self.assertEqual(s['valid_cut_count'],2)
            self.assertEqual(s['valid_proof_count'],1)
            rewritten=(out/s['proofs'][0]['directory']/'model.mps').read_text()
            self.assertIn(' X CUT0 -1\n',rewritten)
            self.assertIn(' Y CUT1 -1\n',rewritten)
            bad=copy.deepcopy(good);bad[-1]['cut_ids'].pop()
            self.replay(Path(td)/'missing-id',model,bad,expected=1)

    def test_unknown_duplicate_and_future_cut_references_fail(self):
        mutations=[]
        e=proof();e['cut_ids']=['MISSING'];mutations.append([root_cut(),e])
        e=proof();e['cut_ids']=['CUT0','CUT0'];mutations.append([root_cut(),e])
        e=proof();e['cut_ids']=[];mutations.append([root_cut(),e])
        mutations.append([proof(),root_cut()])
        e=root_cut();e['certificate']=e['certificate'].replace('PUB_RHS -1','PUB_RHS -2')
        mutations.append([root_cut(),e,proof()])
        mutations.append([root_cut(),root_cut(),proof()])
        mutations.append([root_cut(),proof(),proof()])
        with tempfile.TemporaryDirectory() as td:
            for i,events in enumerate(mutations):
                self.replay(Path(td)/str(i),SCALAR,events,expected=1)

    def test_out_of_box_and_incomplete_box_uses_fail(self):
        events=[]
        for lo,hi in [('-1','3/2'),('5/4','17'),(None,'3/2'),('5/4',None),('2','1')]:
            events.append(proof(lower=lo,upper=hi))
        e=proof();e['bounds']=[];events.append(e)
        e=proof();e['bounds'].append(e['bounds'][0]);events.append(e)
        with tempfile.TemporaryDirectory() as td:
            for i,e in enumerate(events):
                s,_=self.replay(Path(td)/str(i),SCALAR,[root_cut(),e],expected=1)
                self.assertFalse(s['proofs'][0]['independently_valid'])

    def test_tampered_cut_proofs_scopes_coefficients_and_names_fail(self):
        events=[]
        for old,new in [('PUB X -1','PUB X -2'),('PUB_RHS -1','PUB_RHS 0'),
                        ('BETA 3/2','BETA 5/2'),('ROW R U 3','ROW R U 4'),
                        ('PUB X -1','PUB X -1/3'),('BOX X 0 16','BOX X 16 1')]:
            events.append(root_cut(SCALAR_CERT.replace(old,new)))
        for name in ['R','OBJ','bad name','bad\nname','$ignored','']:
            events.append(root_cut(name=name))
        e=root_cut();e['scope']='ORIGINAL_LP';events.append(e)
        e=root_cut();e['coefficient']={'X':'-99'};events.append(e)
        with tempfile.TemporaryDirectory() as td:
            for i,e in enumerate(events):
                s,_=self.replay(Path(td)/str(i),SCALAR,[e,proof()],expected=1)
                self.assertEqual(s['valid_cut_count'],0)

    def test_typed_margins_and_cut_name_laundering_fail_closed(self):
        variants=[]
        e=proof();e.pop('contradiction_margin');e.update(margin_kind='source_enclosure',contradiction_margin_lower_bound='1/8')
        with tempfile.TemporaryDirectory() as td:
            s,_=self.replay(Path(td)/'lower',SCALAR,[root_cut(),e])
            self.assertTrue(s['proofs'][0]['margin_lower_bound_valid'])
            for k,value in [('contradiction_margin','1/2'),('margin_kind','unknown'),('multipliers',{'CUT0':'-1'})]:
                e=proof();e[k]=value;variants.append(e)
            e=proof();e.update(margin_kind='source_enclosure',contradiction_margin_lower_bound='1/8');variants.append(e)
            e=proof();e.pop('contradiction_margin');e.update(margin_kind='source_enclosure',contradiction_margin_lower_bound='1/2');variants.append(e)
            e=proof();e['event']='lp_infeasibility_certificate';variants.append(e)
            e=proof();e['event']='lp_infeasibility_certificate';e.pop('cut_ids');variants.append(e)
            for i,e in enumerate(variants):
                self.replay(Path(td)/str(i),SCALAR,[root_cut(),e],expected=1)

    def test_original_lp_replay_remains_supported_and_unknown_proof_events_fail(self):
        lp=dict(event='lp_infeasibility_certificate',lp_solve=3,
                bounds=[dict(column='X',lower='2',upper='3')],multipliers={'R':'-1'},contradiction_margin='1')
        with tempfile.TemporaryDirectory() as td:
            s,_=self.replay(Path(td)/'lp',SCALAR,[lp])
            self.assertTrue(s['all_lp_certificates_valid'])
            s,_=self.replay(Path(td)/'mixed',SCALAR,[root_cut(),lp])
            self.assertEqual(s['mip_in_box_proof_count'],0)
            for i,event in enumerate([{'event':'unknown_cut_certificate'},{'event':'root_gmi_unknown'},
                                      {'event':{}},{'event':[]},{'event':False}]):
                self.replay(Path(td)/str(i),SCALAR,[root_cut(),event],expected=1)

    def test_duplicate_json_events_empty_policy_and_fresh_output(self):
        with tempfile.TemporaryDirectory() as td:
            self.replay(Path(td)/'empty',SCALAR,[],expected=1)
            self.replay(Path(td)/'allow',SCALAR,[],extra=['--allow-empty'])
            event=json.dumps(root_cut()).replace('"cut_id": "CUT0"','"cut_id": "CUT0", "cut_id": "CUT1"')
            self.replay(Path(td)/'duplicate-key',SCALAR,event,expected=1)
            s,out=self.replay(Path(td)/'run',SCALAR,[root_cut(),proof()])
            before={p.relative_to(out):p.read_bytes() for p in out.rglob('*') if p.is_file()}
            p=subprocess.run([sys.executable,str(REPLAY),'--model',str(out.parent/'original.mps'),
                '--trace',str(out.parent/'trace.jsonl'),'--checker',CHECKER,'--out',str(out)],capture_output=True,text=True)
            self.assertEqual(p.returncode,2)
            self.assertEqual(before,{p.relative_to(out):p.read_bytes() for p in out.rglob('*') if p.is_file()})


if __name__=='__main__':
    result=unittest.main(argv=[sys.argv[0],*rest],exit=False)
    print(json.dumps({'validation_counts':COUNTS},sort_keys=True))
    sys.exit(not result.result.wasSuccessful())
