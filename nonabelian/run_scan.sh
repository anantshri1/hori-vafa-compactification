#!/bin/bash
# LOG 19 GPU batch: window convergence (M, co-moving c) at n=8, then n-collapse at n=8,12,16. x=0.3, Ng=16.
set -x
for c in 0 0.5; do for M in 3 4; do python scan.py 8 0.3 $M 16 $c $c $c; done; done
for n in 12 16; do python scan.py $n 0.3 3 16 0.5 0.5 0.5; python scan.py $n 0.3 3 16 0 0 0; done
