"""Shared helpers for the kp_ (Kepler-Poinsot near-miss) side quest.
Score = (n_pierce, n_faces_crossed, n_crossing_points), lexicographically smaller is better.
Only imports existing modules (k12, check, search, symgroup); never edits them."""
import os, json, itertools, time
from fractions import Fraction as Q
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, HERE)
import k12
import check as checkA


def surface_by_num(num):
    """1-indexed, matches the near_*.npy / LOG.md '#N' convention."""
    L = k12.all59()
    for S in L:
        if S.name.endswith(f'_#{num}'):
            return S
    raise KeyError(num)


def to_rationals(H, grid=10_000_000):
    """Round float plane coeffs to a 1e-7 grid, return as list of 4-tuples of Fraction."""
    out = []
    for row in H:
        out.append(tuple(Q(int(round(float(x) * grid)), grid) for x in row))
    return out


def score(rep):
    """rep: check.py report dict with 'S' (list of (i,c,d)) and 'X' (list of (i,j,k)), already filtered to
    GENUINE piercings only (see filter_genuine_X). Returns (n_pierce, n_faces_crossed, n_crossing_points)."""
    n_pierce = len(rep['X'])
    faces_crossed = sorted({i for (i, c, d) in rep['S']})
    return (n_pierce, len(faces_crossed), len(rep['S']))


def filter_genuine_X(S, H, rep):
    """Per the brief: an X entry (i,j,k) from check.py may be induced only by a self-crossing of face i or j (the test
    point v(i,j,k) then sits ON face k's boundary because edge ki or kj is crossed). Count a GENUINE piercing only when
    v(i,j,k) is not strictly inside edge ki or edge kj. Re-derives points exactly as check.py does."""
    Hq = [tuple(Q(x) for x in h) for h in H]
    n = S.n
    from check import point as _point, strictly_inside_segment as _sis
    P = {}
    for t in itertools.combinations(range(n), 3):
        from check import det3
        if det3(Hq[t[0]][:3], Hq[t[1]][:3], Hq[t[2]][:3]) == 0:
            v = None
        else:
            v = _point(Hq, *t)
        for pm in itertools.permutations(t):
            P[pm] = v

    def edge(i, j):
        a, b = S.elink[tuple(sorted((i, j)))]
        return P[(i, j, a)], P[(i, j, b)]

    def inseg(p, A, B):
        if p is None or p == A or p == B:
            return False
        return _sis(p, A, B)

    genuine = []
    for (i, j, k) in rep['X']:
        p = P[(i, j, k)]
        if inseg(p, *edge(k, i)) or inseg(p, *edge(k, j)):
            continue  # induced only by a face self-crossing at k's boundary -> not genuine
        genuine.append((i, j, k))
    return genuine


def full_report(S, H):
    """Exact-check H (floats or Fractions) against S; returns (ok_embedded, rep, genuine_X, sc)."""
    ok, rep = checkA.check(S, H)
    gX = filter_genuine_X(S, H, rep)
    sc = (len(gX), len(sorted({i for (i, c, d) in rep['S']})), len(rep['S']))
    return ok, rep, gX, sc


def min_edge_ratio(S, H):
    """Minimum edge length / RMS vertex radius from centroid, using exact rationals -> float."""
    Hq = [tuple(Q(x) for x in h) for h in H]
    n = S.n
    P = {}
    for t in itertools.combinations(range(n), 3):
        from check import det3, point as _point
        if det3(Hq[t[0]][:3], Hq[t[1]][:3], Hq[t[2]][:3]) == 0:
            continue
        v = _point(Hq, *t)
        for pm in itertools.permutations(t):
            P[pm] = v
    vtri = {tuple(sorted(t)) for t in S.tris}
    verts = {t: P[t] for t in vtri}
    pts = np.array([[float(x) for x in v] for v in verts.values()])
    centroid = pts.mean(0)
    rms = np.sqrt(((pts - centroid) ** 2).sum(1).mean())
    minlen = 1e18
    for (i, j), (a, b) in S.elink.items():
        A = np.array([float(x) for x in P[(i, j, a)]])
        B = np.array([float(x) for x in P[(i, j, b)]])
        minlen = min(minlen, float(np.linalg.norm(A - B)))
    return minlen / rms, rms


def log_line(msg, path=os.path.join(HERE, 'kp_log.txt')):
    ts = time.strftime('%Y-%m-%d %H:%M:%S')
    line = f'{ts}  {msg}'
    with open(path, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
    print(line)
