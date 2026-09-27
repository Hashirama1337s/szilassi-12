"""Hardness estimate: solve the full dual instance under random cubes of m assumptions on the chi(i,j,k,INF) variables.
usage: python cube_probe.py KEY m trials conflict_budget seed"""
import sys, time, random
import k12, run_sat
from pysat.solvers import Solver
try:
    import psutil; psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception: pass
key, m, trials, budget, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
S = {s.name.split('#')[-1]: s for s in k12.all59()}[key]
bd = run_sat.build('dual', S)
fixed = {abs(c[0]) for c in bd.clauses if len(c) == 1}
pool = [bd.var[b] for b in bd.bases if b[-1] == 12 and bd.var[b] not in fixed]
rng = random.Random(seed)
with Solver(name='cadical195', bootstrap_with=bd.clauses) as s:
    for t in range(trials):
        cube = [v if rng.random() < 0.5 else -v for v in rng.sample(pool, m)]
        s.conf_budget(budget)
        t0 = time.time(); r = s.solve_limited(assumptions=cube)
        print(f'#{key} m={m} trial {t}: {"SAT" if r else ("UNSAT" if r is False else "UNKNOWN")} {time.time() - t0:.1f}s', flush=True)
