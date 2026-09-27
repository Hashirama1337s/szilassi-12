"""Combinatorics of the 59 neighbourly 12-vertex triangulations of the orientable genus-6 surface (Lutz's Manifold Page,
file manifold_lex_d2_n12_o1_g6) and of the 7-vertex torus (the Moebius/Csaszar torus, whose dual is the Szilassi
polyhedron). The dual of such a triangulation T is the polyhedron we want: face i of the dual <-> vertex i of T, dual vertex
<-> triangle of T, dual edge <-> edge of T.  Everything here is 0-indexed."""
import itertools, re, os

HERE = os.path.dirname(os.path.abspath(__file__))


def parse_lutz(path=os.path.join(HERE, 'k12_g6_lutz.txt')):
    txt = open(path, encoding='utf-8').read()
    out = []
    for name, body in re.findall(r'(manifold_\S+?)=\s*(\[\[.*?\]\])', txt, re.S):
        tris = [tuple(int(x) - 1 for x in t) for t in re.findall(r'\[(\d+),(\d+),(\d+)\]', body)]
        out.append((name, tris))
    return out


def torus7():
    """Moebius' 7-vertex torus: triangles {i,i+1,i+3} and {i,i+2,i+3} mod 7 (dual = Szilassi polyhedron)."""
    tris = set()
    for i in range(7):
        tris.add(tuple(sorted((i, (i + 1) % 7, (i + 3) % 7))))
        tris.add(tuple(sorted((i, (i + 2) % 7, (i + 3) % 7))))
    return sorted(tris)


class Surface:
    """A neighbourly triangulated closed orientable surface, with a coherent orientation.
    link[i]  : the neighbours of vertex i in cyclic order (the boundary of dual face i, in order);
    elink[(i,j)] (i<j): the two apexes (a, b) of the triangles on edge ij (the endpoints of dual edge ij);
    otris    : the triangles, coherently oriented."""

    def __init__(self, tris, name=''):
        self.name = name
        self.n = n = 1 + max(max(t) for t in tris)
        self.tris = [tuple(sorted(t)) for t in tris]
        assert len(set(self.tris)) == len(self.tris)
        # every pair of vertices is an edge lying in exactly two triangles (neighbourly closed surface)
        ecount = {}
        for t in self.tris:
            for e in itertools.combinations(t, 2):
                ecount.setdefault(e, []).append(t)
        assert len(ecount) == n * (n - 1) // 2, 'not neighbourly'
        assert all(len(v) == 2 for v in ecount.values()), 'not a closed surface'
        self.elink = {e: tuple(sorted(x for t in ts for x in t if x not in e)) for e, ts in ecount.items()}
        # coherent orientation by breadth-first propagation
        o = {self.tris[0]: self.tris[0]}
        todo = [self.tris[0]]
        while todo:
            t = todo.pop()
            a, b, c = o[t]
            for (x, y) in ((a, b), (b, c), (c, a)):          # oriented edge x->y of t; neighbour must contain y->x
                e = tuple(sorted((x, y)))
                (u,) = [s for s in ecount[e] if s != t]
                z = [w for w in u if w not in e][0]
                want = (y, x, z)
                if u in o:
                    got = o[u]
                    rots = {got, got[1:] + got[:1], got[2:] + got[:2]}
                    assert want in rots, 'non-orientable'
                else:
                    o[u] = want
                    todo.append(u)
        assert len(o) == len(self.tris)
        self.otris = [o[t] for t in self.tris]
        # vertex links as cycles, in the orientation: triangle (i, x, y) contributes the step x -> y
        self.link = {}
        for i in range(n):
            nxt = {}
            for (a, b, c) in self.otris:
                for (p, q, r) in ((a, b, c), (b, c, a), (c, a, b)):
                    if p == i:
                        nxt[q] = r
            cyc = [min(nxt)]
            while nxt[cyc[-1]] != cyc[0]:
                cyc.append(nxt[cyc[-1]])
            assert len(cyc) == n - 1 == len(nxt), 'link is not a single cycle'
            self.link[i] = cyc
        V, E, F = n, len(ecount), len(self.tris)
        self.genus = (2 - (V - E + F)) // 2

    def dual_counts(self):
        return dict(faces=self.n, vertices=len(self.tris), edges=self.n * (self.n - 1) // 2, genus=self.genus)


def all59():
    return [Surface(t, name) for name, t in parse_lutz()]


if __name__ == '__main__':
    S7 = Surface(torus7(), 'torus7')
    print('torus7', S7.dual_counts())
    L = all59()
    print(len(L), 'surfaces;', {tuple(sorted(S.dual_counts().items())) for S in L})
