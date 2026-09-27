"""Validate the OM encodings against geometry: random plane configurations (dual) and random point sets (vertex)."""
import itertools, random, numpy as np, sys
from fractions import Fraction
import k12, check, omsat

def chirotope(vecs):
    chi = {}
    for b in itertools.combinations(range(len(vecs)), 4):
        d = np.linalg.det(np.array([vecs[x] for x in b], dtype=float))
        assert abs(d) > 1e-9
        chi[b] = 1 if d > 0 else -1
    return chi

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
surfs = [k12.Surface(k12.torus7(), 't7')] + k12.all59()[:6]
mism, total, embedded_cases = 0, 0, 0
for trial in range(60):
    S = surfs[trial % len(surfs)]
    n = S.n
    a = rng.normal(size=(n, 3)); a /= np.linalg.norm(a, axis=1)[:, None]
    H = np.hstack([a, (-1 + 0.3 * rng.normal(size=(n, 1)))])
    Hq = [[Fraction(round(x * 10**6), 10**6) for x in h] for h in H]
    Hf = np.array([[float(x) for x in h] for h in Hq])
    chi = chirotope(list(Hf) + [np.array([0, 0, 0, 1.0])])
    assert omsat.gp3_eval(chi, list(range(n + 1))) == 0
    okA, rep = check.check(S, Hq)
    bd = omsat.Builder(n + 1, 'eval', chi)
    v = omsat.dual_constraints(bd, S, n)
    sA, sO = set(rep['S']), set(v['S'])
    # S tests must agree exactly; X tests must agree whenever faces i and j have no S violation (else the test point may
    # sit on a boundary, where the checker counts 'inside' and the parity formula is not defined)
    badS = {f for f, _, _ in rep['S']}
    xA = {x for x in rep['X'] if x[0] not in badS and x[1] not in badS}
    xO = {x for x in v['X'] if x[0] not in badS and x[1] not in badS}
    total += 1
    if sA != sO or xA != xO:
        mism += 1
        print('MISMATCH', S.name[-3:], 'S', len(sA ^ sO), 'X', len(xA ^ xO), sorted(xA ^ xO)[:4])
print('dual: trials', total, 'mismatches', mism)

# vertex version vs a float geometric test (segment-triangle crossing)
def seg_tri(A, B, P, Q, R):
    def o(a, b, c, d): return np.sign(np.linalg.det(np.array([b - a, c - a, d - a])))
    if o(P, Q, R, A) == o(P, Q, R, B): return False
    s1, s2, s3 = o(A, B, P, Q), o(A, B, Q, R), o(A, B, R, P)
    return s1 == s2 == s3
mv = 0
for trial in range(40):
    S = surfs[trial % len(surfs)]; n = S.n
    X = rng.normal(size=(n, 3))
    chi = chirotope([np.append(x, 1.0) for x in X])
    geo = sorted(((a, b), t) for a, b in itertools.combinations(range(n), 2) for t in S.tris
                 if a not in t and b not in t and seg_tri(X[a], X[b], *[X[q] for q in t]))
    bd = omsat.Builder(n, 'eval', chi)
    om = sorted(omsat.vertex_constraints(bd, S))
    mv += geo != om
print('vertex: trials 40 mismatches', mv)
