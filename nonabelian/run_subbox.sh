#!/bin/bash
# LOG 48: coherent sub-box sums (stride 1), exact vs ladder-LG, growing boxes, eps = 0.2. Gate first.
set -e; date -u
python smoke_subbox.py
for hw in 4 8 12; do python subbox.py 0.2 $hw ladder,flip; date -u; done
python subbox.py 0.2 6 far; date -u
