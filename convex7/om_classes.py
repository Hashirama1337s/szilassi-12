"""Descriptive: how many chirotope classes satisfy the base ENC1
instance (no convexity constraint)? Classes are taken up to map automorphisms (face relabelling), mirror image (global
sign change) and reorientation. Bokowski-Schewe 2002 reportedly show the Szilassi embedding is 'essentially unique' via
oriented matroids. If the base instance has a single class, every realisation shares Table 3's oriented matroid.
usage: python om_classes.py [cap] [seconds]  (run 'python convex7.py checks' first; needs python-sat)"""
import sys, os, itertools, time, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import core as c
from pysat.solvers import Solver

cap = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
tlim = float(sys.argv[2]) if len(sys.argv) > 2 else 300
nv, cls, C = c.cnf_enc1(0)
bases = list(itertools.combinations(range(8), 4))
sols = []
t0 = time.time()
with Solver(name='cadical195', bootstrap_with=cls) as s:
    while len(sols) < cap and time.time() - t0 < tlim and s.solve():
        m = s.get_model()
        bits = tuple(1 if m[q] > 0 else -1 for q in range(70))
        sols.append(bits)
        s.add_clause([-(q + 1) if bits[q] > 0 else (q + 1) for q in range(70)])
    exhausted = len(sols) < cap and time.time() - t0 < tlim
print('solutions (one per reorientation class):', len(sols), 'exhausted' if exhausted else 'NOT exhausted',
      f'{time.time() - t0:.0f} s', flush=True)

# map automorphisms
T = set(c.S7.tris)
auts = [p for p in itertools.permutations(range(7)) if {tuple(sorted(p[x] for x in t)) for t in T} == T]
par = c.omsat.perm_parity
idx = {B: q for q, B in enumerate(bases)}


import numpy as np
S = np.array(sols, dtype=np.int8)                                   # (n, 70)
key = {bytes((S[k] > 0).astype(np.uint8)): k for k in range(len(S))}
F = np.array([[(-1) ** sum((R >> b) & 1 for b in B if b < 7) for B in bases] for R in range(128)], dtype=np.int8)
fixed = [idx[B] for B in c.omsat.reorientation_fix(c.omsat.Builder(8, 'cnf'), list(range(8)))]
perms = []
for p in auts:
    pinv = [0] * 8
    for i in range(7):
        pinv[p[i]] = i
    pinv[7] = 7
    src, sg = [], []
    for B in bases:
        t = tuple(pinv[b] for b in B)
        src.append(idx[tuple(sorted(t))]); sg.append((-1) ** par(t))
    perms.append((np.array(src), np.array(sg, dtype=np.int8)))
parent = list(range(len(S)))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
missing = 0
for k in range(len(S)):
    for src, sg in perms:
        for mir in (1, -1):
            v = S[k][src] * sg * mir
            W = F * v                                               # all 128 reorientations
            ok = np.where((W[:, fixed] > 0).all(axis=1))[0]
            w = W[ok[0]]
            j = key.get(bytes((w > 0).astype(np.uint8)))
            if j is None:
                missing += 1; continue
            parent[find(j)] = find(k)
classes = {}
for k in range(len(S)):
    classes.setdefault(find(k), []).append(k)
print('automorphisms of the map:', len(auts), '; images not found in the solution set:', missing)
print('classes up to relabelling, mirror and reorientation:', len(classes), '; sizes:', sorted(len(v) for v in classes.values()))
chi3 = {tuple(json.loads(k.replace('(', '[').replace(')', ']'))): v
        for k, v in json.load(open(os.path.join(c.WORK, 'chi_table3.json'))).items()}
v = np.array([chi3[B] for B in bases], dtype=np.int8)
W = F * v; w = W[np.where((W[:, fixed] > 0).all(axis=1))[0][0]]
j = key.get(bytes((w > 0).astype(np.uint8)))
print('Table 3 (perturbed) found among solutions:', j is not None, '; its class size:', len(classes[find(j)]) if j is not None else None)

# 7-element (projective) restriction: drop the plane at infinity
b7 = [q for q, B in enumerate(bases) if 7 not in B]
F7 = F[:, b7]
canon7 = set()
for k in range(len(S)):
    best = None
    for src, sg in perms:
        for mir in (1, -1):
            v = (S[k][src] * sg * mir)[b7]
            W = F7 * v
            keys = [bytes((w > 0).astype(np.uint8)) for w in W]
            m = min(keys)
            best = m if best is None or m < best else best
    canon7.add(best)
print('classes of the 7-plane (projective) restriction, up to relabelling, mirror, reorientation:', len(canon7))
