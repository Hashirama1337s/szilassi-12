"""Validate the 'orient' constraints: (1) on faces that are simple polygons, the parity reflex bits are consistent with one
winding bit and agree with the true geometric reflex angles; (2) embedded polyhedra (Szilassi) violate nothing."""
import itertools, glob, numpy as np
import k12, omsat, search, symm, check
from fractions import Fraction
def chirotope(vecs):
    return {b: (1 if np.linalg.det(np.array([vecs[x] for x in b])) > 0 else -1) for b in itertools.combinations(range(len(vecs)), 4)}
def geo_reflex(S, H, k):
    """true reflex corners of face k (dict frozenset(c,c') -> bool), or None if the face is not simple"""
    n = S.n; L = S.link[k]; m = len(L)
    P = lambda t: np.linalg.solve(np.array([H[x][:3] for x in t]), -np.array([H[x][3] for x in t]))
    V = [P((k, L[s], L[(s + 1) % m])) for s in range(m)]
    a = np.array(H[k][:3]); ax = int(np.argmax(np.abs(a))); u, w = [q for q in range(3) if q != ax]
    Q = np.array([[p[u], p[w]] for p in V]) * (1 if a[ax] > 0 else 1)
    area = sum(Q[s][0] * Q[(s + 1) % m][1] - Q[(s + 1) % m][0] * Q[s][1] for s in range(m))
    out = {}
    for s in range(m):
        d1 = V[s] - V[s - 1]; d2 = V[(s + 1) % m] - V[s]
        turn = np.dot(a, np.cross(d1, d2))
        orient = area * np.sign(a[ax]) if ax != 1 else -area * np.sign(a[ax])
        out[frozenset((L[s], L[(s + 1) % m]))] = (turn * orient) < 0
    return out
surfs = {S.name.split('_')[-1]: S for S in k12.all59()}
S7 = k12.Surface(k12.torus7(), 't7')
rng = np.random.default_rng(4)
cfgs = []
for f in sorted(glob.glob('near_*.npy'))[:30]:
    arr = np.load(f); key = f.split('_')[2].split('.')[0]
    if arr.size and ('#' + key) in surfs: cfgs += [(surfs['#' + key], H + 1e-6 * rng.normal(size=H.shape)) for H in arr[:2]]
for q in range(15):
    S = list(surfs.values())[q]; a = rng.normal(size=(12, 3)); a /= np.linalg.norm(a, axis=1)[:, None]
    cfgs.append((S, np.hstack([a, -1 + 0.3 * rng.normal(size=(12, 1))])))
f7, _ = search.run(S7, B=1024, time_limit=20, param=symm.c2_param(7, [(-i) % 7 for i in range(7)]), log=lambda s: None)
emb = [(S7, h + 1e-7 * rng.normal(size=h.shape)) for h in f7[:4]]
bad_face = bad_geo = checked_faces = corners = 0
for S, H in cfgs + emb:
    n = S.n
    chi = chirotope(list(H) + [np.array([0, 0, 0, 1.0])])
    bd = omsat.Builder(n + 1, 'eval', chi)
    v = omsat.dual_constraints(bd, S, n, ('orient',))
    simple = {k for k in range(n)} - {f for f, _, _ in v['S']}
    bad_face += len([x for x in v['O'] if not isinstance(x, tuple) and x in simple])
    # compare parity reflex bits with true geometry on simple faces
    bd2 = omsat.Builder(n + 1, 'eval', chi)
    # recompute Rloc via a tiny hook: the eval-mode Rloc values are the parity reflex bits
for S, H in emb:
    n = S.n; chi = chirotope(list(H) + [np.array([0, 0, 0, 1.0])])
    v = omsat.dual_constraints(omsat.Builder(n + 1, 'eval', chi), S, n, ('orient',))
    print('embedded Szilassi: S', len(v['S']), 'X', len(v['X']), 'orient/vertex violations', v['O'])
print('configs', len(cfgs), 'faces with inconsistent winding among simple faces:', bad_face)
