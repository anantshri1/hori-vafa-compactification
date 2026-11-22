#!/bin/bash
# LOG 54: gates (abelian windows, ladder sub-box, partition coverage at small x), then the full partition at eps = 0.3.
set -e; date -u
(cd .. && python window_gate.py)
python smoke_subbox.py
python coverage_gate.py
python partition.py 0.3 0.3 0 0 8 12 2.0 -48; date -u
