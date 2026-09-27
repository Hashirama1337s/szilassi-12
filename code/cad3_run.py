import subprocess, sys, time
cnf, cap = sys.argv[1], float(sys.argv[2])
try:
    import psutil; psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception: pass
t0 = time.time()
try:
    r = subprocess.run([r'tools\cadical\cadical.exe', '-q', cnf], capture_output=True, text=True, timeout=cap)
    res = {10: 'SAT', 20: 'UNSAT'}.get(r.returncode, f'rc{r.returncode}')
except subprocess.TimeoutExpired:
    res = 'UNKNOWN(time limit)'
open('sat_cad3.log', 'a').write(f'{time.ctime()} {cnf}: {res} after {time.time() - t0:.0f}s (cadical 3.0.1)\n')
