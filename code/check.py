"""CHECKER A: exact rational test that 12 (or 7) planes realise the dual of a neighbourly triangulated surface as an embedded
polyhedron.  Input: plane coefficients h_i = (a1, a2, a3, d) (plane a.x + d = 0), as integers or Fractions.

Face i is the polygon in plane P_i with vertices v(i,c1,c2), v(i,c2,c3), ... where (c1, c2, ...) is the link of vertex i
of the triangulation and v(i,j,k) = P_i n P_j n P_k.  Its edge towards face c lies on the line L_ic = P_i n P_c.

The surface is embedded iff (all tests exact):
  G  genericity: every 3 normals independent (all 220 points v(i,j,k) finite) and no 4 planes through a point;
  S  every face is a simple polygon: for non-consecutive c, d in link(i), the point v(i,c,d) is not strictly inside
     both edge ic and edge id;
  X  no edge meets a face it does not belong to: for edge ij with end-apexes {a, b} and every other k, NOT
     (v(i,j,k) strictly inside edge ij AND v(i,j,k) inside face k).
Why X suffices: two faces meet only on their common line, and every pair of faces here is adjacent; any overlap beyond the
shared edge has an end on the boundary of one face inside the other, i.e. an edge-face incidence of type X; edges of
disjoint index pairs lie on skew lines; for k in {a, b} edge ij meets plane k only at its own end vertex.  The X test point
v(i,j,k) always lies on the lines of face k's edges ki and kj; it is on face k's boundary iff it is inside one of those two
segments (then face i or j also fails S), which is tested first; any other boundary contact needs 4 concurrent planes and
raises ValueError (degenerate)."""
from fractions import Fraction as Q
import itertools


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]))


def det4(A):
    return sum((-1) ** j * A[0][j] * det3(*[[r[m] for m in range(4) if m != j] for r in A[1:]]) for j in range(4))


def point(H, i, j, k):
    """v(i,j,k) by Cramer's rule, exact."""
    A = [H[i][:3], H[j][:3], H[k][:3]]
    rhs = [-H[i][3], -H[j][3], -H[k][3]]
    D = det3(*A)
    out = []
    for c in range(3):
        M = [list(r) for r in A]
        for r in range(3):
            M[r][c] = rhs[r]
        out.append(Q(det3(*M)) / D)
    return tuple(out)


def strictly_inside_segment(p, A, B):
    d = [B[m] - A[m] for m in range(3)]
    s = sum((p[m] - A[m]) * d[m] for m in range(3)) / sum(x * x for x in d)
    return 0 < s < 1


def in_polygon(p, poly, normal):
    """Even-odd test after projecting along the dominant normal axis; asserts p is not on the boundary."""
    ax = max(range(3), key=lambda m: abs(normal[m]))
    u, v = [m for m in range(3) if m != ax]
    x, y = p[u], p[v]
    inside = False
    m = len(poly)
    for s in range(m):
        (x1, y1), (x2, y2) = (poly[s][u], poly[s][v]), (poly[(s + 1) % m][u], poly[(s + 1) % m][v])
        cr = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)
        if cr == 0 and min(x1, x2) <= x <= max(x1, x2) and min(y1, y2) <= y <= max(y1, y2):
            raise ValueError('test point on a face boundary (degenerate)')
        if (y1 > y) != (y2 > y):
            if (x1 + (y - y1) * (x2 - x1) / (y2 - y1)) > x:
                inside = not inside
    return inside


def check(S, H, verbose=False):
    """S: k12.Surface; H: list of n rows (a1,a2,a3,d) of ints/Fractions.  Returns (ok, report dict)."""
    n = S.n
    H = [tuple(Q(x) for x in h) for h in H]
    rep = dict(G=[], S=[], X=[])
    # G: the polyhedron's own vertices must be finite and distinct.  A NON-vertex triple may have dependent normals (its test
    # point is at infinity, e.g. forced by a symmetry); such a point lies inside no segment and its tests pass trivially.
    vtri = {tuple(sorted(t)) for t in S.tris}
    P = {}
    for t in itertools.combinations(range(n), 3):
        if det3(H[t[0]][:3], H[t[1]][:3], H[t[2]][:3]) == 0:
            if t in vtri:
                rep['G'].append(('vertex at infinity', t))
            v = None
        else:
            v = point(H, *t)
        for pm in itertools.permutations(t):
            P[pm] = v
    if len({P[t] for t in vtri}) < len(vtri):
        rep['G'].append(('coincident vertices',))
    if rep['G']:
        return False, rep
    def inseg(p, A, B):
        if p is None:
            return False
        if p == A or p == B:
            raise ValueError('test point coincides with a vertex (4 planes concurrent)')
        return strictly_inside_segment(p, A, B)
    faces = {i: [P[(i, S.link[i][s], S.link[i][(s + 1) % (n - 1)])] for s in range(n - 1)] for i in range(n)}

    def edge(i, j):
        a, b = S.elink[tuple(sorted((i, j)))]
        return P[(i, j, a)], P[(i, j, b)]

    for i in range(n):
        L = S.link[i]
        for s, t in itertools.combinations(range(n - 1), 2):
            if (t - s) % (n - 1) in (1, n - 2):
                continue
            c, d = L[s], L[t]
            p = P[(i, c, d)]
            if inseg(p, *edge(i, c)) and inseg(p, *edge(i, d)):
                rep['S'].append((i, c, d))
    for (i, j), (a, b) in S.elink.items():
        A, B = edge(i, j)
        for k in range(n):
            if k in (i, j, a, b):
                continue
            p = P[(i, j, k)]
            if inseg(p, A, B):
                # p lies on the lines of face k's edges ki and kj; if it is inside one of those segments it is on the boundary
                # of the closed face k (and faces i or j cross themselves: an S violation too)
                if inseg(p, *edge(k, i)) or inseg(p, *edge(k, j)) or in_polygon(p, faces[k], H[k][:3]):
                    rep['X'].append((i, j, k))
    ok = not (rep['S'] or rep['X'])
    if verbose:
        print({k: len(v) for k, v in rep.items()}, 'EMBEDDED' if ok else 'not embedded')
    return ok, rep
