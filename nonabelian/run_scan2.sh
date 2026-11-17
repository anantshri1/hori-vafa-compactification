#!/bin/bash
# LOG 23: aliasing-safe (Ng=20) polytope map at n=8, x=0.78: L,R in [-9,6], T in [-4,7], m_B = 0,-2,-4,-6,-8.
set -x; date -u
for mB in -4 0 -8 -2 -6; do python diag.py 8 0.78 $mB -9 6 20 -4 7; date -u; done
