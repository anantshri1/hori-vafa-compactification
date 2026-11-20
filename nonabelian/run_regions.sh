#!/bin/bash
# LOG 47: gate on smoke tests, then weight maps by branch pattern at eps = 0.2, 0.1 (u-box +-3 nu*, stride nu*/4).
set -e; date -u
python smoke_asserts.py
python region_scan.py 0.2 4 3; date -u
python region_scan.py 0.1 4 3; date -u
