"""enc2_geom -- exact (rational) geometry of the dual polyhedron of a neighbourly surface, written independently of
the chirotope encoding.  Faces are projected to a coordinate plane (drop the largest normal component: an affine
bijection of the face plane), and everything is decided with exact 2-D orientation tests, segment tests, a
horizontal-ray point-in-polygon test, and explicit clipping of each face against each line L_ij.

embedded(): every face polygon simple (no two non-adjacent edges meet) AND for every pair of faces i<j the set
(F_i n L_ij) n (F_j n L_ij) equals the common edge (computed by clipping both faces against L_ij and testing every
breakpoint and every open piece between breakpoints).

cond(): the exact geometric truth of every SIMP / PIERCE / REFL / LOC condition of enc2_core, decided here with
coordinates (not with the chirotope)."""
from fractions import Fraction as Fr
import itertools


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def det4(M):
    s = 0
    for c in range(4):
        if M[0][c] == 0:
            continue
        minor = [[M[r][k] for k in range(4) if k != c] for r in range(1, 4)]
        s += (-1) ** c * M[0][c] * det3(minor)
    return s


def orient(p, q, r):
    v = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    return (v > 0) - (v < 0)


def on_seg(p, q, r):
    return (orient(p, q, r) == 0 and min(p[0], q[0]) <= r[0] <= max(p[0], q[0])
            and min(p[1], q[1]) <= r[1] <= max(p[1], q[1]))


def seg_inter(p1, p2, q1, q2):
    o1, o2, o3, o4 = orient(p1, p2, q1), orient(p1, p2, q2), orient(q1, q2, p1), orient(q1, q2, p2)
    if o1 * o2 < 0 and o3 * o4 < 0:
        return True
    return on_seg(p1, p2, q1) or on_seg(p1, p2, q2) or on_seg(q1, q2, p1) or on_seg(q1, q2, p2)


def in_closed(pt, edges):
    """pt in the closed polygon (boundary, or odd crossing number of the ray to +x; half-open rule)."""
    for p, q in edges:
        if on_seg(p, q, pt):
            return True
    cnt = 0
    for p, q in edges:
        if (p[1] > pt[1]) != (q[1] > pt[1]):
            xc = p[0] + (pt[1] - p[1]) * (q[0] - p[0]) / (q[1] - p[1])
            if xc > pt[0]:
                cnt += 1
    return cnt % 2 == 1


class Geo:
    def __init__(self, S, H):
        """S: k12.Surface; H: n planes (a1,a2,a3,d) with rational/integer entries (plane a.x + d = 0)."""
        self.S, self.n = S, S.n
        self.H = [tuple(Fr(x) for x in h) for h in H]
        self._pt = {}
        self.ax = [max(range(3), key=lambda k: abs(h[k])) for h in self.H]
        self.edges3 = {}                       # (i, c) -> 3-D endpoints of the edge of face i toward c
        self.edges2 = {}                       # i -> list of (c, p2, q2) in face i's projection
        for i in range(self.n):
            L = S.link[i]
            m = len(L)
            lst = []
            for t in range(m):
                P, Q = self.pt(i, L[t - 1], L[t]), self.pt(i, L[t], L[(t + 1) % m])
                self.edges3[(i, L[t])] = (P, Q)
                lst.append((L[t], self.proj(i, P), self.proj(i, Q)))
            self.edges2[i] = lst
        self.poly = {i: [(p, q) for _, p, q in self.edges2[i]] for i in range(self.n)}

    # -------- basic --------
    def pt(self, i, j, k):
        key = tuple(sorted((i, j, k)))
        if key not in self._pt:
            A = [self.H[x][:3] for x in key]
            r = [-self.H[x][3] for x in key]
            D = det3(A)
            assert D != 0, 'parallel planes'
            X = []
            for c in range(3):
                Ac = [[r[row] if col == c else A[row][col] for col in range(3)] for row in range(3)]
                X.append(det3(Ac) / D)
            self._pt[key] = tuple(X)
        return self._pt[key]

    def proj(self, i, X):
        a = self.ax[i]
        return tuple(X[k] for k in range(3) if k != a)

    def chi_exact(self, bases):
        """exact chirotope signs of the given 4-tuples (index n = plane at infinity)."""
        rows = list(self.H) + [(Fr(0), Fr(0), Fr(0), Fr(1))]
        out = []
        for B in bases:
            d = det4([rows[b] for b in B])
            out.append((d > 0) - (d < 0))
        return out

    # -------- lines --------
    def line_breaks(self, f, A, B):
        """parameters t (point A + t(B-A), A, B in plane f) where the line meets the boundary of face f."""
        A2, B2 = self.proj(f, A), self.proj(f, B)
        D = (B2[0] - A2[0], B2[1] - A2[1])
        DD = D[0] * D[0] + D[1] * D[1]

        def par(r):
            return ((r[0] - A2[0]) * D[0] + (r[1] - A2[1]) * D[1]) / DD
        out = set()
        for p, q in self.poly[f]:
            op, oq = orient(A2, B2, p), orient(A2, B2, q)
            if op == 0 and oq == 0:
                out.add(par(p))
                out.add(par(q))
            elif op == 0:
                out.add(par(p))
            elif oq == 0:
                out.add(par(q))
            elif op * oq < 0:
                e = (q[0] - p[0], q[1] - p[1])
                num = (p[0] - A2[0]) * e[1] - (p[1] - A2[1]) * e[0]
                den = D[0] * e[1] - D[1] * e[0]
                out.add(num / den)
        return sorted(out)

    def member(self, f, A, B, t):
        X = tuple(A[k] + t * (B[k] - A[k]) for k in range(3))
        return in_closed(self.proj(f, X), self.poly[f])

    def ends(self, i, j):
        a, b = self.S.elink[tuple(sorted((i, j)))]
        return a, b, self.pt(i, j, a), self.pt(i, j, b)

    # -------- conditions (geometric truth) --------
    def simp_bad(self, i, s, t):
        p1, p2 = [self.proj(i, X) for X in self.edges3[(i, s)]]
        q1, q2 = [self.proj(i, X) for X in self.edges3[(i, t)]]
        return seg_inter(p1, p2, q1, q2)

    def reflex(self, i, j, k):
        """points of L_ij just beyond vertex ijk (away from edge ij) belong to closed face i."""
        a, b, A, B = self.ends(i, j)
        br = self.line_breaks(i, A, B)
        if k == a:
            c = [t for t in br if t < 0]
            rep = (max(c) + 0) / 2 if c else Fr(-1)
        else:
            assert k == b
            c = [t for t in br if t > 1]
            rep = (min(c) + 1) / 2 if c else Fr(2)
        return self.member(i, A, B, rep)

    def refl_bad(self, i, j, k):
        return self.reflex(i, j, k) and self.reflex(j, i, k)

    def loc_bad(self, i, j, k):
        return self.reflex(i, j, k) != self.reflex(i, k, j)

    def pierce_bad(self, i, c, j):
        x = self.pt(i, j, c)
        P, Q = self.edges3[(i, c)]
        on_e = on_seg(self.proj(i, P), self.proj(i, Q), self.proj(i, x))
        return on_e and in_closed(self.proj(j, x), self.poly[j])

    def cond_bad(self, fam, key):
        if fam == 'SIMP':
            return self.simp_bad(*key)
        if fam == 'REFL':
            tri, i, j = key
            (k,) = [x for x in tri if x not in (i, j)]
            return self.refl_bad(i, j, k)
        if fam == 'LOC':
            tri, i = key
            j, k = [x for x in tri if x != i]
            return self.loc_bad(i, j, k)
        if fam == 'PIERCE':
            return self.pierce_bad(*key)
        raise KeyError(fam)

    # -------- the embedding test (independent of the condition list) --------
    def simple(self, i):
        E = self.edges2[i]
        m = len(E)
        for s, t in itertools.combinations(range(m), 2):
            if (t - s) % m in (1, m - 1):
                continue
            if seg_inter(E[s][1], E[s][2], E[t][1], E[t][2]):
                return False
        return True

    def pair_ok(self, i, j):
        a, b, A, B = self.ends(i, j)
        T = sorted(set(self.line_breaks(i, A, B)) | set(self.line_breaks(j, A, B)) | {Fr(0), Fr(1)})
        reps = list(T) + [(T[k] + T[k + 1]) / 2 for k in range(len(T) - 1)] + [T[0] - 1, T[-1] + 1]
        for t in reps:
            bi, bj = self.member(i, A, B, t), self.member(j, A, B, t)
            if 0 <= t <= 1:
                assert bi and bj, 'common edge not in both faces?'
            elif bi and bj:
                return False
        return True

    def embedded(self, detail=False):
        simp = [self.simple(i) for i in range(self.n)]
        pairs = {(i, j): self.pair_ok(i, j) for i, j in itertools.combinations(range(self.n), 2)}
        emb = all(simp) and all(pairs.values())
        if detail:
            return emb, sum(not s for s in simp), sum(not p for p in pairs.values())
        return emb
