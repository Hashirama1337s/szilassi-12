import subprocess, sys, time
key = sys.argv[1]
try:
    import psutil; psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception: pass
t0 = time.time()
try:
    r = subprocess.run([r'tools\cadical\cadical.exe', '-q', f'cnf/dual_{key}_full.cnf'], capture_output=True, text=True, timeout=21600)
    res = {10: 'SAT', 20: 'UNSAT'}.get(r.returncode, f'UNKNOWN(rc {r.returncode})')
except subprocess.TimeoutExpired:
    res = 'UNKNOWN(time limit)'
open('sat_full.log', 'a').write(f'{time.ctime()} full #{key}: {res} after {time.time() - t0:.0f}s (cadical 3.0.1, cap 6 h)\n')
if res == 'SAT':
    open(f'model_full_{key}.txt', 'w').write(r.stdout)
