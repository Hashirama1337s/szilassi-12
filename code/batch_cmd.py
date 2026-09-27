"""Run a list of commands (one per line in a file) with WORKERS parallel and a hard time limit each. usage: batch_cmd.py FILE SECS WORKERS TAG"""
import sys, time, subprocess, shlex
from concurrent.futures import ThreadPoolExecutor
f, secs, workers, tag = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
cmds = [l.strip() for l in open(f) if l.strip()]
log = open(f'sat_{tag}.log', 'a', encoding='utf-8'); log.write(f'=== {time.ctime()} {len(cmds)} cmds {secs}s x{workers}\n'); log.flush()
def run(c):
    try:
        r = subprocess.run([sys.executable] + shlex.split(c), capture_output=True, text=True, timeout=secs)
        out = r.stdout.strip() or ('ERROR ' + c + ' ' + r.stderr.strip()[-200:])
    except subprocess.TimeoutExpired:
        out = f'UNKNOWN (time limit {secs:.0f}s): {c}'
    log.write(out + '\n'); log.flush(); return out
with ThreadPoolExecutor(workers) as ex:
    for o in ex.map(run, cmds): print(o, flush=True)
log.write('DONE\n'); log.close()
