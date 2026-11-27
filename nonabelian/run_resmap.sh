#!/bin/bash
# LOG 66: residual map at eps = 0.3 (gate: coverage at small x first).
set -e; date -u
python coverage_gate.py
python resmap.py; date -u
