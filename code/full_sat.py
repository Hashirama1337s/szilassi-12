"""Full (asymmetric) dual instance with extra encodings. usage: python full_sat.py KEY [EXTRA]"""
import sys, time, json
import k12, run_sat
from pysat.solvers import Solver
try:
    import psutil; psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception: pass
key = sys.argv[1]; extra = tuple(sys.argv[2].split(',')) if len(sys.argv) > 2 else ()
S = {s.name.split('#')[-1]: s for s in k12.all59()}[key]
t0 = time.time(); bd = run_sat.build('dual', S, None, extra); tb = time.time() - t0
with Solver(name='cadical195', bootstrap_with=bd.clauses) as s:
    r = s.solve(); model = s.get_model() if r else None
if r:
    m = set(x for x in model if x > 0)
    json.dump({','.join(map(str, b)): (1 if bd.var[b] in m else -1) for b in bd.bases}, open(f'chi_full_{key}.json', 'w'))
print(f'full #{key:3s} extra={",".join(extra) or "-"} {"SAT" if r else "UNSAT"} vars {bd.nv} clauses {len(bd.clauses)} build {tb:.1f}s solve {time.time() - t0 - tb:.1f}s', flush=True)
