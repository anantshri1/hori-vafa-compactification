import subprocess, sys, re
out = subprocess.run([sys.executable, 'subbox.py', '0.2', '1', 'ladder'], capture_output=True, text=True).stdout
print(out.strip())
r = float(re.search(r"\|exact-LG\|/sum\|\.\| (\S+)", out).group(1))
print(f"SMOKE ladder sub-box exact vs LG rel {r:.2e} (need < 3e-2)")
if r >= 3e-2: print("SMOKE FAILED"); sys.exit(1)
print("SMOKE PASSED")
