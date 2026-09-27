import itertools, glob, numpy as np, torch
from fractions import Fraction
import k12, check, omsat, search, symm
def chirotope(vecs):
    chi = {}
    for b in itertools.combinations(range(len(vecs)), 4):
        d = np.linalg.det(np.array([vecs[x] for x in b], dtype=float)); assert abs(d) > 1e-12
        chi[b] = 1 if d > 0 else -1
    return chi
surfs = {S.name.split('_')[-1]: S for S in k12.all59()}
S7 = k12.Surface(k12.torus7(), 't7')
cfgs = []
for f in sorted(glob.glob('near_A_*.npy')):
    key = '#' + f.split('_')[2]
    arr = np.load(f)
    cfgs += [(surfs[key], H) for H in arr[:5]]
f7, info = search.run(S7, B=1024, time_limit=30, noise=0.003, log=lambda s: None, seed=11)
cfgs += [(S7, x[2]) for x in info['near'][:10]] + [(S7, h) for h in f7[:3]]
rng = np.random.default_rng(5)
for q in range(20):
    S = [S7] + list(surfs.values())
    S = S[q % 8]; a = rng.normal(size=(S.n, 3)); a /= np.linalg.norm(a, axis=1)[:, None]
    cfgs.append((S, np.hstack([a, -1 + 0.3 * rng.normal(size=(S.n, 1))])))
tot = dict(S=0, Spos=0, X=0, Xpos=0); mism = 0
for S, H in cfgs:
    n = S.n
    H = H + 1e-7 * np.random.default_rng(0).normal(size=H.shape)       # break exact symmetry (points at infinity)
    Hq = [[Fraction(round(x * 10**7), 10**7) for x in h] for h in H]
    Hf = np.array([[float(x) for x in h] for h in Hq])
    chi = chirotope(list(Hf) + [np.array([0, 0, 0, 1.0])])
    try: okA, rep = check.check(S, Hq)
    except ValueError as e: print('degenerate', e); continue
    v = omsat.dual_constraints(omsat.Builder(n + 1, 'eval', chi), S, n)
    P = {}
    for tt in itertools.combinations(range(n), 3):
        v3 = check.point(Hq, *tt)
        for pm in itertools.permutations(tt): P[pm] = v3
    def seg(i, j):
        a, b = S.elink[tuple(sorted((i, j)))]; return P[(i, j, a)], P[(i, j, b)]
    def onb(i, j, k):   # X test point on the boundary of face k (ambiguous: skip)
        p = P[(i, j, k)]; return check.strictly_inside_segment(p, *seg(k, i)) or check.strictly_inside_segment(p, *seg(k, j))
    badS = set()
    nS = sum(1 for i in range(n) for s, t in itertools.combinations(range(n - 1), 2) if (t - s) % (n - 1) not in (1, n - 2))
    xtests = [(i, j, k) for (i, j), (a, b) in S.elink.items() for k in range(n) if k not in (i, j, a, b) and not onb(i, j, k)]
    tot['S'] += nS; tot['Spos'] += len(rep['S']); tot['X'] += len(xtests)
    xs = set(xtests); xA = {x for x in rep['X'] if x in xs}; xO = {x for x in v['X'] if x in xs}
    tot['Xpos'] += len(xA)
    if set(rep['S']) != set(v['S']) or xA != xO:
        mism += 1; print('MISMATCH', S.name[-3:], sorted(xA ^ xO)[:5])
print('configs', len(cfgs), 'mismatches', mism, tot)
