#!/bin/bash
# LOG 24: contain both L,R edges. n=8, x=0.78, Ng=20; L,R in [-14,12], T in [-5,8]; m_B = 0, -4, -8.
set -x; date -u
for mB in -8 0 -4; do python diag.py 8 0.78 $mB -14 12 20 -5 8; date -u; done
