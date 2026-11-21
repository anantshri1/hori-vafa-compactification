import subprocess, sys, re
out = subprocess.run([sys.executable, 'partition.py', '0.3', '0.45', '0.03', '4', '6', '6', '1.0'], capture_output=True, text=True).stdout
line = [l for l in out.splitlines() if l.startswith('COVERAGE')][0]; print(line)
rel = float(re.search(r'rel (\S+)', line).group(1))
print("COVERAGE GATE", "PASSED" if rel < 1e-5 else "FAILED"); sys.exit(0 if rel < 1e-5 else 1)
