#!/bin/bash
# LOG 60: same full partition (eps = 0.3) with the BRANCH-CORRECT model: diagonals ladder-LG if nu>0 else exact, B,-T exact, vertical 1.
set -e; date -u
(cd .. && python window_gate.py)
python coverage_gate.py
python partition.py 0.3 0.3 0 0 8 12 2.0 -48 branch; date -u
