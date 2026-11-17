#!/bin/bash
# x=0.3 m_B=-4 done (hfjob_diag partial log, 16 Nov). Now x=0.78, batched per-sector evaluation.
set -x
python diag.py 8 0.78 -4 -8 3 16
python diag.py 8 0.78 -7 -10 2 16
