"""SAT experiments on single dual instances: python sat_exp.py SURFACE VARIANT SOLVER SECONDS"""
import sys, time, threading, json
import k12, omsat
from pysat.solvers import Solver
try:
    import psutil; psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception: pass
key, variant, solver, secs = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
S = {s.name.split('#')[-1]: s for s in k12.all59()}[key]
n = S.n
bd = omsat.Builder(n + 1, 'cnf')
omsat.gp3(bd, list(range(n + 1)))
if variant == 'Sonly':
    import itertools
    cache = {}
    def inn(i, j, k):
        a, b = S.elink[tuple(sorted((i, j)))]
        return bd.XOR(bd.neg(i, j, k, a), bd.neg(i, j, a, n), bd.neg(i, j, k, b), bd.neg(i, j, b, n))
    for i in range(n):
        L = S.link[i]; m = len(L)
        for s, t in itertools.combinations(range(m), 2):
            if (t - s) % m in (1, m - 1): continue
            bd.forbid(inn(i, L[s], L[t]), inn(i, L[t], L[s]))
else:
    omsat.dual_constraints(bd, S, n)
omsat.reorientation_fix(bd, list(range(n + 1)))
t0 = time.time()
with Solver(name=solver, bootstrap_with=bd.clauses) as s:
    tm = threading.Timer(secs, s.interrupt); tm.start()
    r = s.solve_limited(expect_interrupt=True); tm.cancel()
print(f'#{key} {variant} {solver} vars {bd.nv} clauses {len(bd.clauses)} -> {r} in {time.time() - t0:.0f}s', flush=True)
