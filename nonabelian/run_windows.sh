#!/bin/bash
# LOG 49: gate (abelian windows vs known 2d answer + ladder sub-box) then Gaussian-windowed mirror sums.
set -e; date -u
(cd .. && python window_gate.py)
python smoke_subbox.py
for s in 0.4 0.6; do python window_mirror.py 0.2 $s ladder,Tflip,Bflip,far; date -u; done
