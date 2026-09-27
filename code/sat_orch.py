import subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
secs = float(sys.argv[1]); exps = [e.split(':') for e in sys.argv[2].split(';')]
log = open('sat_exp.log', 'a', encoding='utf-8'); log.write(f'=== {time.ctime()} {secs}s\n'); log.flush()
def run(e):
    key, var, sol = e
    try:
        r = subprocess.run([sys.executable, 'sat_exp.py', key, var, sol, str(secs)], capture_output=True, text=True, timeout=secs + 300)
        out = r.stdout.strip() or r.stderr.strip()[-300:]
    except subprocess.TimeoutExpired:
        out = f'#{key} {var} {sol} -> TIMEOUT (killed at {secs + 300:.0f}s)'
    log.write(out + '\n'); log.flush(); return out
with ThreadPoolExecutor(len(exps)) as ex:
    for o in ex.map(run, exps): print(o, flush=True)
log.write('DONE\n')
