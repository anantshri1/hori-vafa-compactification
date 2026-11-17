#!/bin/bash
# LOG 21: x=0.78 needs wider, asymmetric boxes (T pushed positive). L,R in [lo,hi], T in [tlo,thi].
set -x
python diag.py 8 0.78 -4 -14 6 16 -4 14
