"""enc2_core -- SECOND, independently derived SAT encoding for the embeddability of the dual polyhedron of a
neighbourly triangulated surface T (faces <-> vertices of T), in terms of the uniform rank-4 chirotope chi on the
face planes h_i = (a_i, d_i) plus h_INF = (0,0,0,1).

Basic predicate (derived from the 4-D cross product): for three planes p = {i,j,k} meeting in the finite point
X = P_i n P_j n P_k and a fourth plane l,
        sign h_l(X) = chi(i,j,k,l) * chi(i,j,k,INF)                       ("side of plane l at point ijk").
Along the line L_ij the affine function h_l vanishes exactly at P_i n P_j n P_l, so "point ijl lies between points
ijk and ijk'" <=> side_l(ijk) != side_l(ijk').  Every condition below is built from side() and XOR/AND/OR only.

Conditions (each a set of clauses over chi and Tseitin gate variables):
  GP     3-term Grassmann-Pluecker, all x<y and 4 others.
  SIMP   (i; s,t): the non-adjacent edges i-s and i-t of face i do not intersect.  The lines meet at point ist;
         it lies on edge i-s iff the endpoints of edge i-s are on opposite sides of plane t (and symmetrically).
  PIERCE (i,c; j): the point x = ijc (where edge i-c of face i crosses line L_ij; c not an apex of edge ij) must not
         lie both on edge i-c and in the closed face j.  "x in face j" = x on edge j-c  OR  x on edge j-i  OR  odd
         number of crossings of a ray from x along L_jc (direction h_i > 0, line pushed to the h_c > 0 side: the
         half-open rule, so the edge j-c lying on L_jc is never counted and its two neighbours count iff their far
         endpoint has h_c > 0).
  REFL   (v=ijk; i,j): faces i and j must not both contain the points of L_ij just beyond the vertex v (away from
         edge ij).  "Face i contains them" = parity of the edges i-d crossing L_ij beyond v (ray along L_ij).
  LOC    (v=ijk; i): the two ways (along L_ij, along L_ik) of saying "face i is reflex at v" agree
         (a tautology for every real generic configuration; kept as a valid extra constraint).
For a realizable uniform chirotope, SIMP & PIERCE & REFL  <=>  embedded (see enc2_README in the report)."""
import itertools, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k12


def perm_sign(t):
    s = 1
    for a, b in itertools.combinations(range(len(t)), 2):
        if t[a] > t[b]:
            s = -s
    return s


# the 16 forbidden sign patterns of a 3-term GP relation (terms t1=x0x1, t2=-x2x3, t3=x4x5 all equal)
GP_BAD = [bits for bits in itertools.product((1, -1), repeat=6)
          if bits[0] * bits[1] == -bits[2] * bits[3] == bits[4] * bits[5]]
assert len(GP_BAD) == 16


class Enc:
    """Build the instance for surface S.  sym: optional list of (perm of the n faces, detM) -- symmetric instance:
    chi(g(B)) = detM * chi(B) with the element INF fixed; implemented by substituting chi variables by orbit reps."""

    def __init__(self, S, fams=('GP', 'SIMP', 'REFL', 'PIERCE', 'LOC'), sym=None):
        self.S, self.n = S, S.n
        n = S.n
        self.INF = n
        self.bases = list(itertools.combinations(range(n + 1), 4))
        self.bidx = {B: k + 1 for k, B in enumerate(self.bases)}
        self.nchi = len(self.bases)
        self.nvar = self.nchi
        self.TRUE = self._new()
        self.defs = [[self.TRUE]]
        self.gates = []
        self.cache = {}
        self.conds = []                        # (family, key, [clauses])
        self.nbr = {}
        for i in range(n):
            L = S.link[i]
            m = len(L)
            self.nbr[i] = {L[t]: (L[t - 1], L[(t + 1) % m]) for t in range(m)}
            for c in L:                        # the neighbours of c in link(i) are the apexes of edge ic
                assert set(self.nbr[i][c]) == set(S.elink[tuple(sorted((i, c)))])
        self.rep = {v: v for v in range(1, self.nchi + 1)}   # chi var -> signed representative literal
        self.contradiction = False
        if sym:
            self._symmetrize(sym)
        self.fams = fams
        self._build()
        self._finalize()

    # ---------------- symmetry ----------------
    def _symmetrize(self, sym):
        par = {v: (v, 1) for v in range(1, self.nchi + 1)}

        def find(v):
            s = 1
            while par[v][0] != v:
                w, t = par[v]
                s *= t
                v = w
            return v, s

        for perm, detM in sym:
            g = list(perm) + [self.INF]
            for B in self.bases:
                gB = tuple(g[b] for b in B)
                # chi(gB as ordered tuple) = detM * chi(B);  chi(gB) = perm_sign(gB) * var(sorted gB)
                u, su = find(self.bidx[tuple(sorted(gB))])
                w, sw = find(self.bidx[B])
                rel = perm_sign(gB) * detM            # var(sorted gB) = rel * var(B)
                if u == w:
                    if su != rel * sw:
                        self.contradiction = True
                else:
                    par[u] = (w, rel * sw * su)
        for v in range(1, self.nchi + 1):
            r, s = find(v)
            self.rep[v] = s * r

    # ---------------- gates ----------------
    def _new(self):
        self.nvar += 1
        return self.nvar

    def chi(self, t):
        return perm_sign(t) * self.rep[self.bidx[tuple(sorted(t))]]

    def XOR(self, a, b):
        s = 1
        if a < 0:
            a, s = -a, -s
        if b < 0:
            b, s = -b, -s
        if a == b:
            return -s * self.TRUE
        if a > b:
            a, b = b, a
        key = ('x', a, b)
        x = self.cache.get(key)
        if x is None:
            x = self._new()
            self.cache[key] = x
            self.defs += [[-x, a, b], [-x, -a, -b], [x, -a, b], [x, a, -b]]
            self.gates.append((x, 0, a, b))
        return s * x

    def AND(self, a, b):
        T = self.TRUE
        if a == b:
            return a
        if a == -b or a == -T or b == -T:
            return -T
        if a == T:
            return b
        if b == T:
            return a
        if abs(a) > abs(b):
            a, b = b, a
        key = ('a', a, b)
        x = self.cache.get(key)
        if x is None:
            x = self._new()
            self.cache[key] = x
            self.defs += [[-x, a], [-x, b], [x, -a, -b]]
            self.gates.append((x, 1, a, b))
        return x

    def OR(self, a, b):
        return -self.AND(-a, -b)

    def XORN(self, lits):
        lits = list(lits)
        if not lits:
            return -self.TRUE
        while len(lits) > 1:
            nxt = [self.XOR(lits[k], lits[k + 1]) for k in range(0, len(lits) - 1, 2)]
            if len(lits) % 2:
                nxt.append(lits[-1])
            lits = nxt
        return lits[0]

    # ---------------- geometric predicates ----------------
    def side(self, p, l):
        """literal: h_l > 0 at the point P_p0 n P_p1 n P_p2."""
        p = tuple(sorted(p))
        return -self.XOR(self.chi(p + (l,)), self.chi(p + (self.INF,)))

    def sdiff(self, p, q, l):
        """literal: points p and q (both 3-sets of planes) lie on opposite sides of plane l."""
        return self.XOR(self.side(p, l), self.side(q, l))

    def crosses(self, i, d, j):
        """literal: edge i-d of face i contains the point P_i n P_d n P_j (its endpoints straddle plane j)."""
        a, b = self.nbr[i][d]
        assert j not in (i, d, a, b)
        return self.sdiff((i, d, a), (i, d, b), j)

    def reflex(self, i, j, k):
        """v = ijk a triangle.  literal: points of L_ij just beyond v (on the side away from edge ij) lie in face i
        (even-odd rule): parity of the edges i-d crossing L_ij on that ray."""
        a, b = self.nbr[i][j]
        assert k in (a, b)
        kk = b if k == a else a                 # other endpoint of edge ij is the point ijkk
        terms = []
        for d in self.S.link[i]:
            if d in (j, a, b):
                continue
            onray = self.sdiff((i, j, d), (i, j, kk), k)     # ijd beyond v: opposite side of plane k from ijkk
            terms.append(self.AND(self.crosses(i, d, j), onray))
        return self.XORN(terms)

    def inface_parts(self, j, c, i):
        """x = P_i n P_j n P_c, a point of L_jc and L_ji inside plane j.  Returns literals (on edge j-c, on edge j-i,
        odd ray parity); x is in closed face j iff their OR is true."""
        p, q = self.nbr[j][c]
        on_jc = self.sdiff((j, c, p), (j, c, q), i)
        a, b = self.nbr[j][i]
        on_ji = self.sdiff((i, j, a), (i, j, b), c)
        terms = []
        for d in self.S.link[j]:
            if d in (c, i):
                continue
            if d in (p, q):                    # edge j-d has one endpoint (jdc) on L_jc: counts iff far end has h_c>0
                (ds,) = [e for e in self.nbr[j][d] if e != c]
                cr = self.side((j, d, ds), c)
            else:
                cr = self.crosses(j, d, c)
            onray = self.side((j, c, d), i)
            terms.append(self.AND(cr, onray))
        return on_jc, on_ji, self.XORN(terms)

    # ---------------- conditions ----------------
    def _build(self):
        S, n, fams = self.S, self.n, self.fams
        if 'GP' in fams:
            N = n + 1
            for x, y in itertools.combinations(range(N), 2):
                rest = [e for e in range(N) if e not in (x, y)]
                for a, b, c, d in itertools.combinations(rest, 4):
                    L = [self.chi((x, y, a, b)), self.chi((x, y, c, d)), self.chi((x, y, a, c)),
                         self.chi((x, y, b, d)), self.chi((x, y, a, d)), self.chi((x, y, b, c))]
                    cls = [[-L[k] if bits[k] == 1 else L[k] for k in range(6)] for bits in GP_BAD]
                    self.conds.append(('GP', (x, y, a, b, c, d), cls))
        if 'SIMP' in fams:
            for i in range(n):
                L = S.link[i]
                m = len(L)
                for s_, t_ in itertools.combinations(range(m), 2):
                    if (t_ - s_) % m in (1, m - 1):
                        continue
                    s, t = L[s_], L[t_]
                    self.conds.append(('SIMP', (i, s, t), [[-self.crosses(i, s, t), -self.crosses(i, t, s)]]))
        if 'REFL' in fams or 'LOC' in fams:
            for tri in S.tris:
                if 'REFL' in fams:
                    for i, j in itertools.combinations(tri, 2):
                        (k,) = [x for x in tri if x not in (i, j)]
                        ri, rj = self.reflex(i, j, k), self.reflex(j, i, k)
                        self.conds.append(('REFL', (tri, i, j), [[-ri, -rj]]))
                if 'LOC' in fams:
                    for i in tri:
                        j, k = [x for x in tri if x != i]
                        r1, r2 = self.reflex(i, j, k), self.reflex(i, k, j)
                        self.conds.append(('LOC', (tri, i), [[-r1, r2], [r1, -r2]]))
        if 'PIERCE' in fams:
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    a, b = self.nbr[i][j]
                    for c in S.link[i]:
                        if c in (j, a, b):
                            continue
                        oe = self.crosses(i, c, j)
                        parts = self.inface_parts(j, c, i)
                        self.conds.append(('PIERCE', (i, c, j), [[-oe, -q] for q in parts]))

    def _finalize(self):
        lev = np.zeros(self.nvar + 1, dtype=np.int64)
        for (x, op, a, b) in self.gates:
            lev[x] = 1 + max(lev[abs(a)], lev[abs(b)])
        groups = {}
        for (x, op, a, b) in self.gates:
            groups.setdefault((int(lev[x]), op), []).append((x, a, b))
        self.levels = []
        for key in sorted(groups):
            arr = np.array(groups[key], dtype=np.int64)
            self.levels.append((key[1], arr[:, 0], np.abs(arr[:, 1]), arr[:, 1] < 0, np.abs(arr[:, 2]), arr[:, 2] < 0))
        self._cl_cache = {}

    def _clause_arrays(self, fams):
        key = tuple(sorted(fams))
        if key not in self._cl_cache:
            idx, cl, ci = [], [], []
            for k, (fam, _, cls) in enumerate(self.conds):
                if fam in fams:
                    idx.append(k)
                    for c in cls:
                        cl.append(c + [0] * (6 - len(c)))
                        ci.append(len(idx) - 1)
            cl = np.array(cl, dtype=np.int64)
            self._cl_cache[key] = (np.array(idx), np.abs(cl), cl < 0, cl == 0, np.array(ci))
        return self._cl_cache[key]

    # ---------------- evaluation on concrete chirotopes ----------------
    def evaluate(self, chivals):
        """chivals: bool array (nchi, B), True = chi(sorted basis) = +1.  Returns all variable values (nvar+1, B)."""
        B = chivals.shape[1]
        V = np.zeros((self.nvar + 1, B), dtype=bool)
        V[1:self.nchi + 1] = chivals
        V[self.TRUE] = True
        for op, out, a, na, b, nb in self.levels:
            va = V[a] ^ na[:, None]
            vb = V[b] ^ nb[:, None]
            V[out] = (va ^ vb) if op == 0 else (va & vb)
        return V

    def violated(self, V, fams):
        """Boolean (ncond_in_fams, B): condition violated (some clause false).  Returns (cond indices, matrix)."""
        idx, av, neg, pad, ci = self._clause_arrays(fams)
        lv = V[av] ^ neg[:, :, None]
        lv[pad] = False
        cfalse = ~lv.any(axis=1)
        out = np.zeros((len(idx), V.shape[1]), dtype=np.int64)
        np.add.at(out, ci, cfalse.astype(np.int64))
        return idx, out > 0

    def defs_ok(self, V):
        """all Tseitin definition clauses satisfied by the evaluated assignment (check of the gate encoding)."""
        L = max(len(c) for c in self.defs)
        cl = np.array([c + [0] * (L - len(c)) for c in self.defs], dtype=np.int64)
        lv = V[np.abs(cl)] ^ (cl < 0)[:, :, None]
        lv[cl == 0] = False
        return bool(lv.any(axis=1).all())

    # ---------------- CNF ----------------
    def cnf(self):
        seen, out = set(), []
        for c in self.defs + [c for _, _, cls in self.conds for c in cls]:
            k = tuple(sorted(set(c)))
            if any(-l in k for l in k):
                continue                                   # tautology
            if k not in seen:
                seen.add(k)
                out.append(list(k))
        return out

    def chi_from_H(self, H):
        """float chirotope of planes H (m, n, 4) [batch m] -> bool (nchi, m)."""
        H = np.asarray(H, dtype=np.float64)
        if H.ndim == 2:
            H = H[None]
        m = H.shape[0]
        Hf = np.concatenate([H, np.broadcast_to(np.array([0, 0, 0, 1.0]), (m, 1, 4))], axis=1)
        idx = np.array(self.bases)
        D = np.linalg.det(Hf[:, idx])                    # (m, nchi)
        return (D > 0).T, np.abs(D).min(axis=1)


def surface(key):
    if key == 't7':
        return k12.Surface(k12.torus7(), 't7')
    for S in k12.all59():
        if S.name.endswith(key):
            return S
    raise KeyError(key)


if __name__ == '__main__':
    import time
    for key in ('t7', '#54'):
        t0 = time.time()
        E = Enc(surface(key))
        from collections import Counter
        print(key, 'vars', E.nvar, 'gates', len(E.gates), 'defs', len(E.defs), 'conds', Counter(f for f, _, _ in E.conds),
              'cnf clauses', len(E.cnf()), '%.1fs' % (time.time() - t0))
