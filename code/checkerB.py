#!/usr/bin/env python3
"""
checkerB.py -- independent exact verifier for embedded polyhedral surfaces.

Usage:
    python checkerB.py file1.json [file2.json ...]

INPUT format (JSON):
    {"vertices": [["x", "y", "z"], ...], "faces": [[i0, i1, ...], ...]}
Coordinates are exact rationals written as strings (e.g. "-3/7" or "12"),
parsed with fractions.Fraction. Each face is a cyclic list of vertex
indices; its vertices are supposed to be coplanar and the face is the
(possibly non-convex) simple polygon bounded by that closed cycle.

All decisions are made in exact rational arithmetic. No floating point
is used anywhere in the decision logic.

For every input file this prints one line:
    <case>: EMBEDDED V=.. E=.. F=.. genus=..
    <case>: NOT EMBEDDED V=.. E=.. F=.. genus=.. | faces (i,j) reason: ...
    <case>: DEGENERATE [V=.. E=.. F=.. genus=..] reason: ...
"""

import sys
import os
import json
from fractions import Fraction as F
from collections import defaultdict

# ----------------------------------------------------------------------
# basic exact vector algebra (3D)
# ----------------------------------------------------------------------

def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scal(a, t):
    return (a[0] * t, a[1] * t, a[2] * t)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def is_zero_vec(v):
    return v[0] == 0 and v[1] == 0 and v[2] == 0


def sign(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


def parse_point(p):
    return (F(p[0]), F(p[1]), F(p[2]))


# ----------------------------------------------------------------------
# exact 2D predicates (used after projecting a planar face to 2D)
# ----------------------------------------------------------------------

def orient(a, b, c):
    """sign of the cross product (b-a) x (c-a) in 2D."""
    return sign((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def on_segment_2d(a, b, p):
    """assumes p is collinear with a,b; checks it's within the bounding box."""
    return (min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and
            min(a[1], b[1]) <= p[1] <= max(a[1], b[1]))


def segments_intersect_2d(p1, p2, p3, p4):
    """True if closed segments p1p2 and p3p4 share any point."""
    d1 = orient(p3, p4, p1)
    d2 = orient(p3, p4, p2)
    d3 = orient(p1, p2, p3)
    d4 = orient(p1, p2, p4)
    if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and \
       ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)):
        return True
    if d1 == 0 and on_segment_2d(p3, p4, p1):
        return True
    if d2 == 0 and on_segment_2d(p3, p4, p2):
        return True
    if d3 == 0 and on_segment_2d(p1, p2, p3):
        return True
    if d4 == 0 and on_segment_2d(p1, p2, p4):
        return True
    return False


def point_in_closed_triangle_2d(a, b, c, p):
    d1 = orient(a, b, p)
    d2 = orient(b, c, p)
    d3 = orient(c, a, p)
    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (has_neg and has_pos)


def cross2(ux, uy, vx, vy):
    return ux * vy - uy * vx


# ----------------------------------------------------------------------
# ear-clipping triangulation of a simple (possibly non-convex) polygon
# ----------------------------------------------------------------------

class EarClipError(Exception):
    pass


def triangulate(face_indices, proj2d):
    """
    face_indices: list of ORIGINAL vertex indices in cyclic order (len n>=3)
    proj2d: dict original_index -> (x,y) Fraction tuple, the projected
            2D coordinates (guaranteed non-degenerate projection of a
            planar simple polygon).
    Returns a list of triangles, each a tuple of 3 ORIGINAL vertex indices.
    """
    remaining = list(face_indices)
    n0 = len(remaining)
    area2 = F(0)
    for i in range(n0):
        x1, y1 = proj2d[remaining[i]]
        x2, y2 = proj2d[remaining[(i + 1) % n0]]
        area2 += x1 * y2 - x2 * y1
    orient_sign = sign(area2)
    if orient_sign == 0:
        raise EarClipError("zero-area polygon")

    triangles = []
    guard = 0
    max_guard = n0 * n0 * 4 + 10
    while len(remaining) > 3:
        guard += 1
        if guard > max_guard:
            raise EarClipError("ear clipping failed to converge")
        n = len(remaining)
        found = False
        for i in range(n):
            ip = remaining[(i - 1) % n]
            ic = remaining[i]
            inx = remaining[(i + 1) % n]
            a, b, c = proj2d[ip], proj2d[ic], proj2d[inx]
            o = orient(a, b, c)
            if o != orient_sign:
                continue  # reflex or degenerate corner, not an ear
            ear_ok = True
            for k in remaining:
                if k == ip or k == ic or k == inx:
                    continue
                if point_in_closed_triangle_2d(a, b, c, proj2d[k]):
                    ear_ok = False
                    break
            if ear_ok:
                triangles.append((ip, ic, inx))
                remaining.pop(i)
                found = True
                break
        if not found:
            raise EarClipError("no ear found (polygon may not be simple)")
    triangles.append(tuple(remaining))
    return triangles


# ----------------------------------------------------------------------
# exact triangle-triangle intersection in 3D
# ----------------------------------------------------------------------

def triangle_normal(pts):
    return cross(sub(pts[1], pts[0]), sub(pts[2], pts[0]))


def signed_dist(n, p0, p):
    return dot(n, sub(p, p0))


def clip_triangle_against_plane(tri_pts, n, p0):
    """
    Returns ('empty', None) / ('coplanar', None) / ('point', P) /
    ('segment', (A,B)) describing tri_pts ∩ plane(n,p0), exactly.
    """
    d = [signed_dist(n, p0, p) for p in tri_pts]
    s = [sign(x) for x in d]
    pos = [i for i in range(3) if s[i] > 0]
    neg = [i for i in range(3) if s[i] < 0]
    zero = [i for i in range(3) if s[i] == 0]

    if len(zero) == 3:
        return ('coplanar', None)
    if len(pos) == 3 or len(neg) == 3:
        return ('empty', None)
    if len(zero) == 2:
        i, j = zero
        return ('segment', (tri_pts[i], tri_pts[j]))
    if len(zero) == 1:
        i = zero[0]
        others = [k for k in range(3) if k != i]
        if s[others[0]] == s[others[1]]:
            return ('point', tri_pts[i])
        a, b = others
        da, db = d[a], d[b]
        t = da / (da - db)
        cross_pt = add(tri_pts[a], scal(sub(tri_pts[b], tri_pts[a]), t))
        return ('segment', (tri_pts[i], cross_pt))
    # len(zero) == 0, split 1-2
    if len(pos) == 1:
        lone, others = pos[0], neg
    else:
        lone, others = neg[0], pos
    a, b = others
    dl = d[lone]
    da, db = d[a], d[b]
    t1 = dl / (dl - da)
    pt1 = add(tri_pts[lone], scal(sub(tri_pts[a], tri_pts[lone]), t1))
    t2 = dl / (dl - db)
    pt2 = add(tri_pts[lone], scal(sub(tri_pts[b], tri_pts[lone]), t2))
    return ('segment', (pt1, pt2))


def project_2d(p3, drop_axis):
    axes = [0, 1, 2]
    axes.remove(drop_axis)
    return (p3[axes[0]], p3[axes[1]])


def unproject_2d(pt2d, drop_axis, n, p0):
    axes = [0, 1, 2]
    axes.remove(drop_axis)
    coords = [None, None, None]
    coords[axes[0]] = pt2d[0]
    coords[axes[1]] = pt2d[1]
    s = sum(n[a] * (coords[a] - p0[a]) for a in axes)
    coords[drop_axis] = p0[drop_axis] - s / n[drop_axis]
    return tuple(coords)


def sh_clip_convex(subject, clip):
    """Sutherland-Hodgman clip of convex polygon `subject` by convex
    polygon `clip`, both lists of exact 2D points. Returns list of 2D
    points (possibly empty / degenerate)."""

    def area2(poly):
        s = F(0)
        n = len(poly)
        for i in range(n):
            x1, y1 = poly[i]
            x2, y2 = poly[(i + 1) % n]
            s += x1 * y2 - x2 * y1
        return s

    if area2(clip) < 0:
        clip = clip[::-1]

    output = list(subject)
    m = len(clip)
    for i in range(m):
        if not output:
            break
        c1 = clip[i]
        c2 = clip[(i + 1) % m]
        inp = output
        output = []
        k = len(inp)
        for j in range(k):
            cur = inp[j]
            prev = inp[j - 1]
            cur_in = orient(c1, c2, cur) >= 0
            prev_in = orient(c1, c2, prev) >= 0
            if cur_in:
                if not prev_in:
                    output.append(seg_intersect_line(prev, cur, c1, c2))
                output.append(cur)
            else:
                if prev_in:
                    output.append(seg_intersect_line(prev, cur, c1, c2))
    return output


def seg_intersect_line(p, q, a, b):
    abx, aby = b[0] - a[0], b[1] - a[1]
    apx, apy = a[0] - p[0], a[1] - p[1]
    aqpx, aqpy = q[0] - p[0], q[1] - p[1]
    denom = cross2(abx, aby, aqpx, aqpy)
    numer = cross2(abx, aby, apx, apy)
    t = numer / denom
    x = p[0] + t * (q[0] - p[0])
    y = p[1] + t * (q[1] - p[1])
    return (x, y)


def classify_2d_region(poly):
    """poly: list of exact 2D points (output of sh_clip_convex).
    Returns ('empty',None) / ('point',P) / ('segment',(A,B)) / ('polygon',pts)."""
    seen = []
    for p in poly:
        if not seen or seen[-1] != p:
            seen.append(p)
    if len(seen) >= 2 and seen[0] == seen[-1]:
        seen.pop()
    pts = []
    for p in seen:
        if p not in pts:
            pts.append(p)
    if len(pts) == 0:
        return ('empty', None)
    if len(pts) == 1:
        return ('point', pts[0])
    area2 = F(0)
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        area2 += x1 * y2 - x2 * y1
    if area2 != 0:
        return ('polygon', pts)
    x0, y0 = pts[0]
    vx, vy = pts[1][0] - x0, pts[1][1] - y0
    if vx == 0 and vy == 0:
        # pts[1] coincides with pts[0]; find another distinct point for direction
        vx = vy = None
        for p in pts[1:]:
            if p != pts[0]:
                vx, vy = p[0] - x0, p[1] - y0
                break
        if vx is None:
            return ('point', pts[0])

    def keyfun(p):
        return (p[0] - x0) * vx + (p[1] - y0) * vy

    pts_sorted = sorted(pts, key=keyfun)
    a, b = pts_sorted[0], pts_sorted[-1]
    if a == b:
        return ('point', a)
    return ('segment', (a, b))


def triangle_triangle_intersection(tri1, tri2):
    """
    tri1, tri2: lists of 3 exact 3D points.
    Returns ('empty', None) / ('point', P) / ('segment', (A,B)) /
            ('polygon', reason_str)
    """
    n1 = triangle_normal(tri1)
    n2 = triangle_normal(tri2)
    p0_1, p0_2 = tri1[0], tri2[0]

    cr = cross(n1, n2)
    if is_zero_vec(cr):
        # parallel planes (or same plane)
        if dot(n1, sub(p0_2, p0_1)) != 0:
            return ('empty', None)  # parallel, distinct planes
        # coplanar case: exact 2D convex-polygon intersection
        absn = [abs(n1[0]), abs(n1[1]), abs(n1[2])]
        drop = absn.index(max(absn))
        t1_2d = [project_2d(p, drop) for p in tri1]
        t2_2d = [project_2d(p, drop) for p in tri2]
        clipped = sh_clip_convex(t1_2d, t2_2d)
        kind, val = classify_2d_region(clipped)
        if kind == 'empty':
            return ('empty', None)
        if kind == 'point':
            return ('point', unproject_2d(val, drop, n1, p0_1))
        if kind == 'segment':
            a2, b2 = val
            return ('segment', (unproject_2d(a2, drop, n1, p0_1),
                                 unproject_2d(b2, drop, n1, p0_1)))
        if kind == 'polygon':
            return ('polygon', 'coplanar faces overlap over a 2D region')
        raise RuntimeError("unreachable")

    # planes distinct and not parallel: intersect in a line L
    chord1 = clip_triangle_against_plane(tri1, n2, p0_2)
    chord2 = clip_triangle_against_plane(tri2, n1, p0_1)

    if chord1[0] == 'empty' or chord2[0] == 'empty':
        return ('empty', None)
    if chord1[0] == 'coplanar' or chord2[0] == 'coplanar':
        # shouldn't happen for non-degenerate triangles with distinct planes
        raise RuntimeError("degenerate triangle encountered in tri-tri test")

    dirL = cr  # direction of L = plane1 ∩ plane2
    D = dot(dirL, dirL)

    # Reference point for the parametrization MUST lie on L itself.
    # chord1 is triangle1 ∩ plane2, and triangle1 ⊆ plane1, so chord1 ⊆
    # plane1 ∩ plane2 = L exactly -- its first point is a valid ref point.
    # (tri1[0]/tri2[0] are NOT generally on L, unlike a prior buggy draft.)
    refP = chord1[1] if chord1[0] == 'point' else chord1[1][0]

    def as_interval(chord):
        if chord[0] == 'point':
            P = chord[1]
            t = dot(sub(P, refP), dirL)
            return [t, t]
        else:
            A, B = chord[1]
            tA = dot(sub(A, refP), dirL)
            tB = dot(sub(B, refP), dirL)
            return sorted([tA, tB])

    lo1, hi1 = as_interval(chord1)
    lo2, hi2 = as_interval(chord2)
    lo = max(lo1, lo2)
    hi = min(hi1, hi2)
    if lo > hi:
        return ('empty', None)

    def point_at(t):
        s = t / D
        return add(refP, scal(dirL, s))

    if lo == hi:
        return ('point', point_at(lo))
    return ('segment', (point_at(lo), point_at(hi)))


# ----------------------------------------------------------------------
# allowed-region membership (shared vertex / shared edge / nothing)
# ----------------------------------------------------------------------

def point_on_segment_3d(A, B, P):
    AB = sub(B, A)
    AP = sub(P, A)
    if not is_zero_vec(cross(AB, AP)):
        return False
    denom = dot(AB, AB)
    if denom == 0:
        return A == P
    t = dot(AP, AB)
    return 0 <= t <= denom


def point_allowed(P, shared_verts, shared_edges, positions):
    for v in shared_verts:
        if positions[v] == P:
            return True
    for (u, v) in shared_edges:
        if point_on_segment_3d(positions[u], positions[v], P):
            return True
    return False


def segment_allowed(A, B, shared_edges, positions):
    for (u, v) in shared_edges:
        U, V = positions[u], positions[v]
        if point_on_segment_3d(U, V, A) and point_on_segment_3d(U, V, B):
            return True
    return False


def intersection_allowed(kind, val, shared_verts, shared_edges, positions):
    if kind == 'empty':
        return True
    if kind == 'point':
        return point_allowed(val, shared_verts, shared_edges, positions)
    if kind == 'segment':
        A, B = val
        if A == B:
            return point_allowed(A, shared_verts, shared_edges, positions)
        return segment_allowed(A, B, shared_edges, positions)
    if kind == 'polygon':
        return False
    return False


# ----------------------------------------------------------------------
# vertex link (single-cycle) check
# ----------------------------------------------------------------------

def vertex_link_ok(v, faces):
    wedges = []
    for face in faces:
        n = len(face)
        try:
            i = face.index(v)
        except ValueError:
            continue
        prev_v = face[(i - 1) % n]
        next_v = face[(i + 1) % n]
        wedges.append((prev_v, next_v))
    if not wedges:
        return False, "vertex is isolated (not used by any face)"
    deg = defaultdict(int)
    adj = defaultdict(list)
    for (p, q) in wedges:
        deg[p] += 1
        deg[q] += 1
        adj[p].append(q)
        adj[q].append(p)
    nodes = set(deg.keys())
    bad = [n for n in nodes if deg[n] != 2]
    if bad:
        return False, f"neighbor {bad[0]} has link-degree {deg[bad[0]]} (expected 2)"
    start = next(iter(nodes))
    visited = {start}
    stack = [start]
    while stack:
        cur = stack.pop()
        for nb in adj[cur]:
            if nb not in visited:
                visited.add(nb)
                stack.append(nb)
    if visited != nodes:
        return False, "link forms more than one cycle (disconnected)"
    return True, len(wedges)


# ----------------------------------------------------------------------
# main per-case analysis
# ----------------------------------------------------------------------

class CaseResult:
    def __init__(self):
        self.status = None       # 'EMBEDDED' / 'NOT EMBEDDED' / 'DEGENERATE'
        self.V = None
        self.E = None
        self.F = None
        self.genus = None
        self.reason = None
        self.pair = None


def analyze(data):
    res = CaseResult()
    verts_raw = data['vertices']
    faces = data['faces']
    positions = [parse_point(p) for p in verts_raw]
    V = len(positions)
    F_ = len(faces)
    res.V, res.F = V, F_

    # ---- combinatorial pass (independent of geometry) ----
    edge_faces = defaultdict(list)
    for fid, face in enumerate(faces):
        n = len(face)
        if n < 2:
            continue
        for i in range(n):
            u, v = face[i], face[(i + 1) % n]
            edge_faces[frozenset((u, v))].append(fid)
    E = len(edge_faces)
    res.E = E

    manifold_reason = None
    for e, fl in edge_faces.items():
        if len(fl) != 2:
            u, v = tuple(e) if len(e) == 2 else (list(e)[0], list(e)[0])
            manifold_reason = f"edge ({u},{v}) belongs to {len(fl)} face(s) (expected 2)"
            break

    unused = None
    used = set()
    for face in faces:
        used.update(face)
    missing = set(range(V)) - used
    if missing:
        unused = min(missing)
        if manifold_reason is None:
            manifold_reason = f"vertex {unused} is not used by any face"

    if manifold_reason is None:
        for v in range(V):
            ok_link, info = vertex_link_ok(v, faces)
            if not ok_link:
                manifold_reason = f"vertex {v} link is invalid: {info}"
                break

    connected = None
    if manifold_reason is None:
        adjv = defaultdict(set)
        for e in edge_faces:
            if len(e) == 2:
                u, v = tuple(e)
                adjv[u].add(v)
                adjv[v].add(u)
        if V > 0:
            seenv = {0}
            stack = [0]
            while stack:
                cur = stack.pop()
                for nb in adjv[cur]:
                    if nb not in seenv:
                        seenv.add(nb)
                        stack.append(nb)
            connected = (len(seenv) == V)
        else:
            connected = True
        if not connected:
            manifold_reason = f"surface is disconnected ({len(seenv)}/{V} vertices reached from vertex 0)"

    chi = V - E + F_
    if connected and chi % 2 == 0:
        res.genus = (2 - chi) // 2
    else:
        res.genus = f"undefined (V-E+F={chi})"

    # ---- per-face geometric validity + triangulation ----
    face_tris = [None] * F_
    geom_reason = None
    for fid, face in enumerate(faces):
        n = len(face)
        if n < 3:
            geom_reason = f"face {fid} has fewer than 3 vertices"
            break
        if len(set(face)) != n:
            geom_reason = f"face {fid} has a repeated vertex"
            break
        pts = [positions[i] for i in face]

        normal = None
        base = pts[0]
        for a in range(1, n - 1):
            for b in range(a + 1, n):
                e1 = sub(pts[a], base)
                e2 = sub(pts[b], base)
                cr = cross(e1, e2)
                if not is_zero_vec(cr):
                    normal = cr
                    break
            if normal is not None:
                break
        if normal is None:
            geom_reason = f"face {fid} vertices are all collinear"
            break

        offplane = None
        for k in range(n):
            if dot(normal, sub(pts[k], base)) != 0:
                offplane = face[k]
                break
        if offplane is not None:
            geom_reason = f"face {fid} is not planar (vertex {offplane} off-plane)"
            break

        collinear_bad = None
        for i in range(n):
            a = pts[(i - 1) % n]
            b = pts[i]
            c = pts[(i + 1) % n]
            if is_zero_vec(cross(sub(b, a), sub(c, b))):
                collinear_bad = face[i]
                break
        if collinear_bad is not None:
            geom_reason = f"face {fid} has collinear consecutive edges at vertex {collinear_bad}"
            break

        absn = [abs(normal[0]), abs(normal[1]), abs(normal[2])]
        drop = absn.index(max(absn))
        axes = [0, 1, 2]
        axes.remove(drop)
        proj = {face[i]: (pts[i][axes[0]], pts[i][axes[1]]) for i in range(n)}

        self_int = None
        for i in range(n):
            a1, a2 = face[i], face[(i + 1) % n]
            for j in range(i + 1, n):
                if j == i + 1 or (i == 0 and j == n - 1):
                    continue
                b1, b2 = face[j], face[(j + 1) % n]
                p1, p2 = proj[a1], proj[a2]
                p3, p4 = proj[b1], proj[b2]
                if segments_intersect_2d(p1, p2, p3, p4):
                    self_int = (a1, a2, b1, b2)
                    break
            if self_int:
                break
        if self_int:
            geom_reason = (f"face {fid} is self-intersecting "
                            f"(edges ({self_int[0]},{self_int[1]}) and "
                            f"({self_int[2]},{self_int[3]}))")
            break

        try:
            face_tris[fid] = triangulate(face, proj)
        except EarClipError as ex:
            geom_reason = f"face {fid} triangulation failed: {ex}"
            break

    if geom_reason is not None:
        res.status = 'DEGENERATE'
        res.reason = geom_reason
        return res

    if manifold_reason is not None:
        res.status = 'DEGENERATE'
        res.reason = manifold_reason
        return res

    # ---- embedding check: exact triangle-triangle intersection ----
    violation = None
    for f1 in range(F_):
        face1 = faces[f1]
        n1 = len(face1)
        set1 = set(face1)
        edges1 = set(frozenset((face1[i], face1[(i + 1) % n1])) for i in range(n1))
        for f2 in range(f1 + 1, F_):
            face2 = faces[f2]
            n2 = len(face2)
            set2 = set(face2)
            edges2 = set(frozenset((face2[i], face2[(i + 1) % n2])) for i in range(n2))
            shared_verts = set1 & set2
            shared_edge_sets = edges1 & edges2
            shared_edges = [tuple(e) for e in shared_edge_sets]

            found_here = None
            for t1 in face_tris[f1]:
                tri1_pts = [positions[i] for i in t1]
                for t2 in face_tris[f2]:
                    tri2_pts = [positions[i] for i in t2]
                    kind, val = triangle_triangle_intersection(tri1_pts, tri2_pts)
                    if kind == 'empty':
                        continue
                    if intersection_allowed(kind, val, shared_verts, shared_edges, positions):
                        continue
                    if kind == 'polygon':
                        reason = val
                    elif kind == 'segment':
                        reason = ("triangles intersect along a segment not contained "
                                  f"in a shared edge (triangle {t1} of face {f1} vs "
                                  f"triangle {t2} of face {f2})")
                    else:  # point
                        reason = ("triangles touch at a point that is not a shared "
                                  f"vertex (triangle {t1} of face {f1} vs "
                                  f"triangle {t2} of face {f2})")
                    found_here = reason
                    break
                if found_here:
                    break
            if found_here:
                violation = (f1, f2, found_here)
                break
        if violation:
            break

    if violation:
        f1, f2, reason = violation
        res.status = 'NOT EMBEDDED'
        res.pair = (f1, f2)
        res.reason = reason
        return res

    res.status = 'EMBEDDED'
    return res


def format_result(name, res):
    genus_str = res.genus if isinstance(res.genus, str) else str(res.genus)
    base = f"{name}: {res.status}"
    if res.V is not None:
        base += f" V={res.V} E={res.E} F={res.F} genus={genus_str}"
    if res.status == 'NOT EMBEDDED':
        base += f" | faces {res.pair} reason: {res.reason}"
    elif res.status == 'DEGENERATE':
        base += f" | reason: {res.reason}"
    return base


def main():
    if len(sys.argv) < 2:
        print("usage: python checkerB.py file1.json [file2.json ...]", file=sys.stderr)
        sys.exit(1)
    for path in sys.argv[1:]:
        name = os.path.splitext(os.path.basename(path))[0]
        try:
            with open(path, 'r') as f:
                data = json.load(f)
            res = analyze(data)
        except Exception as ex:
            res = CaseResult()
            res.status = 'DEGENERATE'
            res.reason = f"error while processing: {ex!r}"
        print(format_result(name, res))


if __name__ == '__main__':
    main()
