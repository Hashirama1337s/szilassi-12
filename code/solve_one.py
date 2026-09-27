"""Solve one OM SAT instance (no time limit here; the caller kills it). usage: python solve_one.py MODE KEY"""
import sys, time, json
import k12, run_sat
from pysat.solvers import Solver
try:
    import psutil; psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception: pass
mode, key = sys.argv[1], sys.argv[2]
S = k12.Surface(k12.torus7(), 't7') if key == 't7' else {s.name.split('#')[-1]: s for s in k12.all59()}[key]
t0 = time.time(); bd = run_sat.build(mode, S); tb = time.time() - t0
with Solver(name='cadical195', bootstrap_with=bd.clauses) as s:
    r = s.solve(); model = s.get_model() if r else None; st = s.accum_stats()
if r:
    m = set(x for x in model if x > 0)
    json.dump({','.join(map(str, b)): (1 if bd.var[b] in m else -1) for b in bd.bases}, open(f'chi_{mode}_{key}.json', 'w'))
print(f'{mode:6s} #{key:3s} {"SAT" if r else "UNSAT":7s} vars {bd.nv:7d} clauses {len(bd.clauses):8d} build {tb:5.1f}s solve {time.time() - t0 - tb:8.1f}s conflicts {st.get("conflicts", 0)}', flush=True)
