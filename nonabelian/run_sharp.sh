#!/bin/bash
# LOG 25: edge sharpening. Same observable at (n,x) = (8,0.78) and (12,0.85); m_B = 0 and -n.
set -x; date -u
python diag2d.py 8 0.78 0 -16 12 20 0,2,4
python diag2d.py 8 0.78 -8 -16 12 20 0,2,4
python diag2d.py 12 0.85 0 -22 18 24 0,2,4
python diag2d.py 12 0.85 -12 -22 18 24 0,2,4
date -u
