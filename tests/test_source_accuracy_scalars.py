#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys,tempfile
from fractions import Fraction
from test_source_accuracy import case_text
assert Fraction('1.000000000000000001')-1==Fraction(1,10**18)
assert (Fraction('1.000000000000000001')-1)*10**15==Fraction(1,1000)
with tempfile.TemporaryDirectory(prefix='source-scalar-controls-') as temporary:
 p=Path(temporary);coupled=p/'coupled.mps';row=p/'row.mps';bound=p/'bound.mps';unknown=p/'unknown.mps'
 coupled.write_text(case_text('coupled',large=True,fixed=True)[0].replace('OBJSENSE','OBJNAME OBJ\nOBJSENSE').replace('RHS\n','RHS\n rhs OBJ -0.0001\n'));row.write_text(case_text('row',large=True)[0])
 base=case_text('objective',large=True)[0]
 bound.write_text(base.replace(' X OBJ 1.000000000000000001\n X OBJ -1',' X OBJ 0').replace(' LO b X 0',' LO b X 1.000000000000000001'))
 unknown.write_text(base.replace(' X OBJ 1.000000000000000001\n X OBJ -1',' X OBJ 1.'+'0'*1024+'1'))
 subprocess.run([sys.argv[1],str(coupled),str(row),str(bound),str(unknown)],check=True)
