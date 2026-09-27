"""Validate 'ival': with the chirotope FIXED to a real configuration, the CNF must be satisfiable exactly when the
polyhedron is embedded (checker A). Compared against the base encoding on the same configurations."""
import itertools, glob, numpy as np
from fractions import Fraction
from pysat.solvers import Solver
import k12, omsat, check, search, symm
def chirotope(vecs):
    return {b: (1 if np.linalg.det(np.array([vecs[x] for x in b])) > 0 else -1) for b in itertools.combinations(range(len(vecs)), 4)}
def status(S, chi, extra):
    n = S.n; bd = omsat.Builder(n + 1, 'cnf'); omsat.gp3(bd, list(range(n + 1)))
    omsat.dual_constraints(bd, S, n, extra)
    units = [[bd.var[b] if s > 0 else -bd.var[b]] for b, s in chi.items()]
    with Solver(name='cadical195', bootstrap_with=bd.clauses + units) as s:
        return s.solve()
surfs = {S.name.split('_')[-1]: S for S in k12.all59()}; S7 = k12.Surface(k12.torus7(), 't7')
rng = np.random.default_rng(21); cfgs = []
f1, _ = search.run(S7, B=1024, time_limit=15, param=symm.c2_param(7, [(-i) % 7 for i in range(7)]), log=lambda s: None, seed=3)
cfgs += [(S7, h + 1e-7 * rng.normal(size=h.shape)) for h in f1[:4]]
f2, info = search.run(S7, B=1024, time_limit=20, noise=0.003, log=lambda s: None, seed=4)
cfgs += [(S7, x[2]) for x in info['near'][:6]]
for f in sorted(glob.glob('near_*.npy'))[:12]:
    arr = np.load(f); key = '#' + f.split('_')[2].split('.')[0]
    if arr.size and key in surfs: cfgs.append((surfs[key], arr[0] + 1e-6 * rng.normal(size=arr[0].shape)))
agree = mism = emb = 0
for S, H in cfgs:
    n = S.n; chi = chirotope(list(H) + [np.array([0, 0, 0, 1.0])])
    Hq = [[Fraction(round(x * 10**9), 10**9) for x in h] for h in H]
    try: okA = check.check(S, Hq)[0]
    except ValueError: continue
    b0 = status(S, chi, ()); b1 = status(S, chi, ('orient',)); b2 = status(S, chi, ('orient', 'ival'))
    emb += okA
    if b0 == b1 == b2 == okA: agree += 1
    else: mism += 1; print('MISMATCH', S.name[-4:], 'checkerA', okA, 'base', b0, 'orient', b1, 'ival', b2)
print('configs', agree + mism, 'embedded', emb, 'agree', agree, 'mismatch', mism)
