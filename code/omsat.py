"""Oriented-matroid SAT encodings.

DUAL (face) version, the open problem: planes h_0..h_{n-1} in R^4 (plane a.x + d = 0) plus h_inf = (0,0,0,1); the rank-4
chirotope chi on these n+1 elements (uniform, WLOG: the realisation space is open).  With chi(i,j,k,INF) = sign det(a_i,a_j,a_k):
  order on line L_ij (direction a_i x a_j):  point l after point k   <=>  chi(ijlk) chi(ijl INF) chi(ijk INF) = +
  v(i,j,k) strictly inside edge ij (apexes a, b)    <=>  chi(ijka) chi(ija INF) chi(ijkb) chi(ijb INF) = -
  side of v(k,a,c) w.r.t. plane i                  =    chi(kaci) chi(kac INF)
Constraints (the exact counterparts of checker A's tests S and X):
  S(i,c,d): not( in(ic,d) and in(id,c) )                    for c, d non-consecutive in link(i)
  X(ij,k):  not( in(ij,k) and inside(v(i,j,k), face k) )    for k not in {i, j, a, b}
  inside(p = v(k,i,j), face k): parity of boundary crossings on the ray from p along L_ki in direction a_k x a_i:
     transversal crossings at v(k,i,c) with in(kc,i), c not in {k,i,j,a',b'}, c after j;
     plus the collinear edge ki itself (apexes a', b') iff it lies after p and the polygon leaves it to opposite sides.
VERTEX version (Schewe's control): points p_0..p_{n-1} (homogeneous (p,1)), acyclic uniform chirotope; for every edge ab and
every triangle pqr disjoint from it: not( chi(pqra) != chi(pqrb) and chi(abpq) = chi(abqr) = chi(abrp) ).

Every predicate is built once as an expression over sign bits, so the SAME code both writes CNF (Tseitin) and evaluates on a
concrete chirotope (used to validate the encoding against checker A on random geometric configurations)."""
import itertools


def perm_parity(t):
    t = list(t); p = 0
    for x in range(len(t)):
        for y in range(x + 1, len(t)):
            if t[x] > t[y]:
                p ^= 1
    return p


class Builder:
    """Boolean circuit builder. Mode 'cnf': literals are DIMACS ints; mode 'eval': literals are Python bools."""

    def __init__(self, N, mode='cnf', chi=None):
        self.N, self.mode, self.chi = N, mode, chi
        self.bases = list(itertools.combinations(range(N), 4))
        self.var = {b: q + 1 for q, b in enumerate(self.bases)}
        self.nv = len(self.bases)
        self.clauses = []
        self.marks = []
        self.T = self.new() if mode == 'cnf' else True
        if mode == 'cnf':
            self.clauses.append([self.T])

    def new(self):
        self.nv += 1
        return self.nv

    # sign bit: TRUE means chi(tuple) = -   (a 'negativity' bit, so products are XORs)
    def neg(self, *t):
        if len(set(t)) < 4:
            raise ValueError(t)
        s = tuple(sorted(t)); par = perm_parity(t)
        if self.mode == 'eval':
            return (self.chi[s] < 0) ^ bool(par)
        v = self.var[s]
        return -v if par == 0 else v            # chi(sorted)=+ <-> var true; neg bit = not var, flipped by parity

    def NOT(self, a):
        return (not a) if self.mode == 'eval' else -a

    def XOR(self, *xs):
        xs = list(xs)
        if self.mode == 'eval':
            r = False
            for x in xs:
                r ^= x
            return r
        if not xs:
            return -self.T
        acc = xs[0]
        for x in xs[1:]:
            y = self.new()
            self.clauses += [[-y, acc, x], [-y, -acc, -x], [y, -acc, x], [y, acc, -x]]
            acc = y
        return acc

    def AND(self, *xs):
        if self.mode == 'eval':
            return all(xs)
        y = self.new()
        for x in xs:
            self.clauses.append([-y, x])
        self.clauses.append([y] + [-x for x in xs])
        return y

    def forbid(self, *xs):
        """assert not (x1 and x2 and ...); in eval mode return True iff the constraint is VIOLATED."""
        if self.mode == 'eval':
            return all(xs)
        self.clauses.append([-x for x in xs])


# ---------------------------------------------------------------------------------------------------------------------
def gp3(bd, elems):
    """3-term Grassmann-Pluecker relations for a uniform rank-4 chirotope on elems (cnf mode only)."""
    for x, y in itertools.combinations(elems, 2):
        rest = [e for e in elems if e not in (x, y)]
        for a, b, c, d in itertools.combinations(rest, 4):
            # sigma1 = chi(xyab)chi(xycd), sigma2 = -chi(xyac)chi(xybd), sigma3 = chi(xyad)chi(xybc): not all equal.
            v = [bd.neg(x, y, a, b), bd.neg(x, y, c, d), bd.neg(x, y, a, c), bd.neg(x, y, b, d), bd.neg(x, y, a, d),
                 bd.neg(x, y, b, c)]
            # negativity of sigma1 = v0^v1 ; sigma2 = 1^v2^v3 ; sigma3 = v4^v5. Forbid all three equal: enumerate the 16
            # assignments of the 6 bits that make them equal and add one blocking clause for each.
            for bits in itertools.product((0, 1), repeat=6):
                s1 = bits[0] ^ bits[1]; s2 = 1 ^ bits[2] ^ bits[3]; s3 = bits[4] ^ bits[5]
                if s1 == s2 == s3:
                    bd.clauses.append([(-v[q] if bits[q] else v[q]) for q in range(6)])


def gp3_eval(chi, elems):
    bad = 0
    for x, y in itertools.combinations(elems, 2):
        rest = [e for e in elems if e not in (x, y)]
        for a, b, c, d in itertools.combinations(rest, 4):
            def s(*t):
                return chi[tuple(sorted(t))] * (-1) ** perm_parity(t)
            s1 = s(x, y, a, b) * s(x, y, c, d); s2 = -s(x, y, a, c) * s(x, y, b, d); s3 = s(x, y, a, d) * s(x, y, b, c)
            bad += (s1 == s2 == s3)
    return bad


def dual_constraints(bd, S, INF, extra=()):
    """Adds (cnf) or evaluates (eval) the S and X constraints of the dual problem. Returns lists of violated tests in eval."""
    n = S.n
    elink = lambda i, j: S.elink[tuple(sorted((i, j)))]
    cache = {}

    def after(i, j, k, l):          # on line ij: l after k
        key = ('af', i, j, k, l)
        if key not in cache:
            cache[key] = bd.NOT(bd.XOR(bd.neg(i, j, l, k), bd.neg(i, j, l, INF), bd.neg(i, j, k, INF)))
        return cache[key]

    def inn(i, j, k):               # v(i,j,k) strictly inside edge ij
        key = ('in', i, j, k)
        if key not in cache:
            a, b = elink(i, j)
            cache[key] = bd.XOR(bd.neg(i, j, k, a), bd.neg(i, j, a, INF), bd.neg(i, j, k, b), bd.neg(i, j, b, INF))
        return cache[key]

    def side_neg(k, a, c, i):       # negativity of the side of v(k,a,c) w.r.t. plane i
        return bd.XOR(bd.neg(k, a, c, i), bd.neg(k, a, c, INF))

    def other_nb(k, a, not_this):   # the neighbour of a in link(k) other than not_this
        L = S.link[k]; m = len(L); q = L.index(a)
        u, w = L[(q - 1) % m], L[(q + 1) % m]
        return w if u == not_this else u

    def inside(k, i, j):            # v(k,i,j) inside face k, by the ray along L_ki
        a, b = elink(k, i)
        terms = []
        for c in range(n):
            if c in (k, i, j, a, b):
                continue
            terms.append(bd.AND(inn(k, c, i), after(k, i, j, c)))
        aa, bb = other_nb(k, a, i), other_nb(k, b, i)
        cross = bd.XOR(side_neg(k, a, aa, i), side_neg(k, b, bb, i))
        terms.append(bd.AND(cross, after(k, i, j, a)))
        return bd.XOR(*terms)

    extra_viol = dict(T=[], R=[])
    if 'trans' in extra:            # explicit transitivity of the order on every line (implied by GP; a shortcut)
        for i, j in itertools.combinations(range(n), 2):
            pts = [c for c in range(n) if c not in (i, j)]
            for k, l, m in itertools.permutations(pts, 3):
                if bd.forbid(after(i, j, k, l), after(i, j, l, m), bd.NOT(after(i, j, k, m))):
                    extra_viol['T'].append((i, j, k, l, m))
    if 'ray2' in extra:             # insideness of v(k,i,j) computed along L_ki and along L_kj must agree (Jordan),
        for k in range(n):          # whenever the point is not on face k's boundary
            for i, j in itertools.combinations([x for x in range(n) if x != k], 2):
                if k in S.elink[tuple(sorted((i, j)))]:
                    continue        # v(k,i,j) is a vertex of face k
                A, B = inside(k, i, j), inside(k, j, i)
                onb1, onb2 = inn(k, i, j), inn(k, j, i)
                if bd.mode == 'cnf':
                    bd.clauses += [[onb1, onb2, -A, B], [onb1, onb2, A, -B]]
                elif not onb1 and not onb2 and A != B:
                    extra_viol['R'].append((k, i, j))
    if 'orient' in extra:
        # Face orientation W_k (negativity of the polygon's winding w.r.t. a_k) + local turns T at each corner give the
        # reflex bit locally: R = T xor W.  Consistency with the global (parity) reflex test, and at most one reflex
        # corner per polyhedron vertex (two reflex corners at a degree-3 vertex make the faces overlap).
        Rloc, orient_bad, Wvar = {}, [], {}
        for k in range(n):
            L = S.link[k]; m = len(L)
            W = bd.new() if bd.mode == 'cnf' else None
            vals = []
            Wvar[k] = W
            for s in range(m):
                c, cp = L[s], L[(s + 1) % m]
                cprev, cnext = L[(s - 1) % m], L[(s + 2) % m]
                sig_in = after(k, c, cprev, cp)          # travel along L_kc (towards v(k,c,cp)) is in +u_kc direction
                sig_out = after(k, cp, c, cnext)         # travel along L_kcp (away from v) is in +u_kcp direction
                T = bd.XOR(bd.NOT(sig_in), bd.NOT(sig_out), bd.neg(k, c, cp, INF))   # negativity of the turn
                terms = []
                for x in range(n):
                    if x in (k, c, cp, cprev):
                        continue
                    beyond = bd.NOT(bd.XOR(after(k, c, cp, x), sig_in))
                    terms.append(bd.AND(inn(k, x, c), beyond))
                Rpar = bd.XOR(*terms)                    # points just past v along L_kc are inside face k
                if bd.mode == 'cnf':
                    R = bd.XOR(T, W)
                    bd.clauses += [[-R, Rpar], [R, -Rpar]]
                    Rloc[(k, frozenset((c, cp)))] = R
                else:
                    vals.append(T != Rpar)               # = W, must be constant along the face
                    Rloc[(k, frozenset((c, cp)))] = Rpar
            if bd.mode != 'cnf':
                Wvar[k] = vals[0] if vals else False
                if len(set(vals)) > 1:
                    orient_bad.append(k)
        for tri in S.tris:
            corners = [Rloc[(f, frozenset(x for x in tri if x != f))] for f in tri]
            for p, q in itertools.combinations(corners, 2):
                if bd.forbid(p, q):
                    orient_bad.append(('V', tri))
        extra_viol['O'] = orient_bad
        bd.Rloc = Rloc                                   # per-corner reflex literals, used by convex7/ (no CNF change)
    if 'jordan3' in extra:
        # A line parallel to L_ij, pushed slightly into the quadrant (si, sj) (sides of P_i, P_j), misses faces i and j
        # and meets face k only near v(i,j,k); a closed embedded surface is crossed an even number of times.
        jbad = []
        def side_neg(k, a, c, i):   # negativity of the side of v(k,a,c) w.r.t. plane i
            return bd.XOR(bd.neg(k, a, c, i), bd.neg(k, a, c, INF))
        def travel(k, i):           # negativity of face k's travel direction along edge ki (w.r.t. u_ki = a_k x a_i)
            L = S.link[k]; m = len(L); q = L.index(i)
            return bd.NOT(after(k, i, L[(q - 1) % m], L[(q + 1) % m]))
        for i, j in itertools.combinations(range(n), 2):
            for si in (False, True):            # negativity bits of the quadrant
                for sj in (False, True):
                    cross = []
                    for k in range(n):
                        if k in (i, j):
                            continue
                        ea = S.elink[tuple(sorted((i, j)))]
                        if k in ea:              # v(i,j,k) is a corner of face k
                            jp = [x for x in S.elink[tuple(sorted((k, i)))] if x != j][0]
                            ip = [x for x in S.elink[tuple(sorted((k, j)))] if x != i][0]
                            s_j = side_neg(k, i, jp, j)          # far end of edge ki: side of P_j
                            s_i = side_neg(k, j, ip, i)          # far end of edge kj: side of P_i
                            inquad = bd.AND(bd.NOT(bd.XOR(s_i, si if bd.mode == 'eval' else (bd.T if si else -bd.T))),
                                            bd.NOT(bd.XOR(s_j, sj if bd.mode == 'eval' else (bd.T if sj else -bd.T))))
                            R = Rloc[(k, frozenset((i, j)))]
                            cross.append(bd.XOR(inquad, R))
                        else:
                            onki, onkj = inn(k, i, j), inn(k, j, i)
                            # interior side of edge ki (as a negativity bit of the side of P_i) = -travel * winding
                            int_i = bd.NOT(bd.XOR(travel(k, i), Wvar[k]))
                            int_j = bd.NOT(bd.XOR(travel(k, j), Wvar[k]))
                            ci = bd.NOT(bd.XOR(int_i, si if bd.mode == 'eval' else (bd.T if si else -bd.T)))
                            cj = bd.NOT(bd.XOR(int_j, sj if bd.mode == 'eval' else (bd.T if sj else -bd.T)))
                            inner = inside(k, i, j)
                            # exactly one of: on edge ki, on edge kj, neither (S forbids both)
                            cross.append(bd.XOR(bd.AND(onki, ci), bd.AND(onkj, cj),
                                                bd.AND(bd.NOT(onki), bd.NOT(onkj), inner)))
                    par = bd.XOR(*cross)
                    if bd.mode == 'cnf':
                        bd.clauses.append([-par])
                    elif par:
                        jbad.append((i, j, si, sj))
        extra_viol['J'] = jbad
    IVAL = {}
    if 'ival' in extra:
        # Inside-status of face f along each of its lines L_fg, propagated between CONSECUTIVE points instead of by
        # parity: IN-(x)/IN+(x) = status just before/after point v(f,g,x) (direction u_fg).  A transversal crossing of
        # edge fx flips it to 'entering' (local, from the winding bit); the collinear edge fg is 'on' the boundary and
        # the status just outside its ends is the reflex bit of those corners; far ends are outside.
        assert 'orient' in extra and bd.mode == 'cnf'
        TT = bd.T
        def eq(x, y, guard):             # guard -> (x == y)
            bd.clauses += [guard + [-x, y], guard + [x, -y]]
        def travel(k, i):
            L = S.link[k]; m = len(L); q = L.index(i)
            return bd.NOT(after(k, i, L[(q - 1) % m], L[(q + 1) % m]))
        for f in range(n):
            for g in range(n):
                if g == f:
                    continue
                pts = [x for x in range(n) if x not in (f, g)]
                a, b = S.elink[tuple(sorted((f, g)))]
                INm = {x: bd.new() for x in pts}
                INp = {x: bd.new() for x in pts}
                IVAL[(f, g)] = INm
                in_edge = {}
                for x in pts:
                    if x in (a, b):
                        in_edge[x] = TT
                        continue
                    on = inn(f, g, x)                        # v(f,g,x) strictly inside edge fg
                    in_edge[x] = on
                    cr = inn(f, x, g)                        # edge fx crosses L_fg at v(f,g,x)
                    interior = bd.NOT(bd.XOR(travel(f, x), Wvar[f]))      # interior side of edge fx (side of P_x)
                    enter = bd.NOT(bd.XOR(bd.neg(f, g, x, INF), interior))
                    eq(INp[x], enter, [-cr])
                    eq(INm[x], -enter, [-cr])
                    eq(INp[x], INm[x], [cr, on])
                Ra, Rb = Rloc[(f, frozenset((g, a)))], Rloc[(f, frozenset((g, b)))]
                a_first = after(f, g, a, b)
                eq(INm[a], Ra, [-a_first]); eq(INp[b], Rb, [-a_first])
                eq(INm[b], Rb, [a_first]); eq(INp[a], Ra, [a_first])
                for p in pts:
                    bd.clauses.append([after(f, g, r, p) for r in pts if r != p] + [-INm[p]])   # first point
                    bd.clauses.append([after(f, g, p, r) for r in pts if r != p] + [-INp[p]])   # last point
                    for q in pts:
                        if q == p:
                            continue
                        btw = []
                        for r in pts:
                            if r in (p, q):
                                continue
                            z = bd.new()
                            bd.clauses += [[-z, after(f, g, p, r)], [-z, after(f, g, r, q)]]
                            btw.append(z)
                        guard = [-after(f, g, p, q)] + btw
                        # consecutive p < q, not both on the collinear edge: status carries over
                        both = bd.AND(in_edge[p], in_edge[q])
                        eq(INm[q], INp[p], guard + [both])
    viol = dict(S=[], X=[], **extra_viol)
    for i in range(n):
        L = S.link[i]; m = len(L)
        for s, t in itertools.combinations(range(m), 2):
            if (t - s) % m in (1, m - 1):
                continue
            c, d = L[s], L[t]
            bd.marks.append((len(bd.clauses), ('S', i, c, d)))
            if bd.forbid(inn(i, c, d), inn(i, d, c)):
                viol['S'].append((i, c, d))
    for (i, j), (a, b) in S.elink.items():
        for k in range(n):
            if k in (i, j, a, b):
                continue
            bd.marks.append((len(bd.clauses), ('X', i, j, k)))
            if bd.forbid(inn(i, j, k), IVAL[(k, i)][j] if 'ival' in extra else inside(k, i, j)):
                viol['X'].append((i, j, k))
    bd.marks.append((len(bd.clauses), ('END',)))
    return viol


def vertex_constraints(bd, S):
    viol = []
    for a, b in itertools.combinations(range(S.n), 2):
        for (p, q, r) in S.tris:
            if a in (p, q, r) or b in (p, q, r):
                continue
            opp = bd.XOR(bd.neg(p, q, r, a), bd.neg(p, q, r, b))
            e1 = bd.NOT(bd.XOR(bd.neg(a, b, p, q), bd.neg(a, b, q, r)))
            e2 = bd.NOT(bd.XOR(bd.neg(a, b, q, r), bd.neg(a, b, r, p)))
            if bd.forbid(opp, e1, e2):
                viol.append(((a, b), (p, q, r)))
    return viol


def acyclic(bd, elems):
    """no positive circuit: for each 5-set, the signs (-1)^m chi(set minus x_m) are not all equal (cnf mode)."""
    for five in itertools.combinations(elems, 5):
        lits = []
        for m in range(5):
            rest = [e for q, e in enumerate(five) if q != m]
            nb = bd.neg(*rest)
            lits.append(-nb if m % 2 else nb)       # negativity of (-1)^m chi(rest)... m odd flips
        # forbid all negativities equal (all True or all False)
        bd.clauses.append([-x for x in lits]); bd.clauses.append(list(lits))


def reorientation_fix(bd, elems):
    """Reorienting a set R of elements multiplies chi(B) by (-1)^|B cap R|; the kernel is {0, all} (|B| = 4 is even), so the
    reorientation classes form GF(2)^E / <ones> of dimension |E|-1.  Fixing chi = + on |E|-1 bases with GF(2)-independent
    incidence vectors picks exactly one chirotope per class (the constraints are reorientation invariant)."""
    E = list(elems); idx = {e: q for q, e in enumerate(E)}
    piv, picked = {}, []
    for bset in itertools.combinations(E, 4):
        w = sum(1 << idx[e] for e in bset)
        while w:
            p = w.bit_length() - 1
            if p in piv:
                w ^= piv[p]
            else:
                piv[p] = w; picked.append(bset); break
        if len(picked) == len(E) - 1:
            break
    for bset in picked:
        bd.clauses.append([bd.var[tuple(sorted(bset))]])
    return picked