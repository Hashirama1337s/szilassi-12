"""Export a plane configuration as a plain polyhedron (exact rational vertices + faces as cyclic vertex lists)."""
import json, itertools
from fractions import Fraction as Q
import check

def polyhedron(S, H):
    H = [tuple(Q(x) for x in h) for h in H]
    vid = {}
    verts = []
    for t in S.tris:
        vid[tuple(sorted(t))] = len(verts)
        verts.append(check.point(H, *t))
    faces = [[vid[tuple(sorted((i, S.link[i][s], S.link[i][(s + 1) % (S.n - 1)])))] for s in range(S.n - 1)] for i in range(S.n)]
    return dict(vertices=[[str(x) for x in v] for v in verts], faces=faces)
