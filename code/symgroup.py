"""Symmetric parametrisations. For a subgroup G of Aut(T) acting freely on the faces, and a faithful representation
rho: G -> O(3) (orientation-preserving automorphisms -> rotations, orientation-reversing -> improper isometries), a
G-invariant configuration is fixed by one plane per face orbit: plane(g r) = (rho(g) a_r, d_r).  Mirror images realise the
same (unoriented) complex, so one representation per abstract group suffices."""
import json, os, itertools, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AUT = json.load(open(os.path.join(HERE, 'automorphisms.json')))

def compose(p, q):            # (p o q)(x) = p[q[x]]
    return tuple(p[x] for x in q)

def order(p):
    e = tuple(range(len(p))); q, k = p, 1
    while q != e:
        q = compose(p, q); k += 1
    return k

def rot(axis, ang):
    x, y, z = np.asarray(axis, float) / np.linalg.norm(axis)
    c, s = math.cos(ang), math.sin(ang)
    return np.array([[c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
                     [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
                     [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)]])

def generate(gens, n):
    """BFS closure; returns dict element -> word (list of generator indices)."""
    e = tuple(range(n)); G = {e: []}; frontier = [e]
    while frontier:
        new = []
        for g in frontier:
            for k, s in enumerate(gens):
                h = compose(s, g)
                if h not in G:
                    G[h] = [k] + G[g]; new.append(h)
        frontier = new
    return G

def representation(gens, mats, n):
    G = generate(gens, n)
    rho = {}
    for g, w in G.items():
        M = np.eye(3)
        for k in reversed(w):
            M = mats[k] @ M
        rho[g] = M
    for g in G:                                          # verify homomorphism on the whole group
        for h in G:
            if not np.allclose(rho[g] @ rho[h], rho[compose(g, h)], atol=1e-9):
                return None
    return rho

def elements(name):
    return [(tuple(p), pres, o) for p, pres, o in AUT[name]['elems']]

S4MAT = np.array([[0.0, -1, 0], [1, 0, 0], [0, 0, -1]])   # rotoreflection: 90 deg about z, then z -> -z

def candidate_groups(name, n=12):
    """Yield (label, rho) for the useful symmetry groups of surface `name`."""
    E = elements(name)
    out, seen = [], set()
    def add(label, gens, mats):
        rho = representation(gens, mats, n)
        key = frozenset(rho) if rho else None
        if rho and key not in seen and all(sum(g[v] == v for v in range(n)) == 0 for g in rho if g != tuple(range(n))):
            seen.add(key); out.append((label, rho))
    inv = [p for p, pres, o in E if o == 2 and pres]
    o3 = [p for p, pres, o in E if o == 3 and pres]
    o4r = [p for p, pres, o in E if o == 4 and not pres]
    # tetrahedral T = A4: s (order 2), t (order 3) with s t of order 3
    for s in inv:
        for t in o3:
            if order(compose(s, t)) == 3 and len(generate([s, t], n)) == 12:
                add('T', [s, t], [np.diag([1.0, -1, -1]), np.array([[0.0, 0, 1], [1, 0, 0], [0, 1, 0]])])
    # dihedral D3: s (order 2), t (order 3), s t of order 2
    for s in inv:
        for t in o3:
            if order(compose(s, t)) == 2 and len(generate([s, t], n)) == 6:
                add('D3', [s, t], [np.diag([1.0, -1, -1]), rot([0, 0, 1], 2 * math.pi / 3)])
    # D2 = V4
    for a, b in itertools.combinations(inv, 2):
        if compose(a, b) == compose(b, a) and len(generate([a, b], n)) == 4:
            add('D2', [a, b], [np.diag([1.0, -1, -1]), np.diag([-1.0, 1, -1])])
    for g in o4r:
        add('S4', [g], [S4MAT])
    for t in o3:
        add('C3', [t], [rot([0, 0, 1], 2 * math.pi / 3)])
    for s in inv:
        add('C2', [s], [np.diag([-1.0, -1, 1])])
    return out

def param(rho, n=12, dev='cuda'):
    """(init_fn, map_fn, reps) for the run() search: Z has shape (B, m, 4), one plane per face orbit."""
    import torch
    reps, face_src = [], {}
    for f in range(n):
        if f in face_src:
            continue
        reps.append(f)
        for g, M in rho.items():
            face_src[g[f]] = (len(reps) - 1, M)
    m = len(reps)
    q = torch.tensor([face_src[f][0] for f in range(n)], device=dev)
    M = torch.tensor(np.stack([face_src[f][1] for f in range(n)]), dtype=torch.float32, device=dev)   # (n,3,3)
    def map_fn(Z):
        a = torch.einsum('fij,bfj->bfi', M, Z[:, q, :3])
        return torch.cat([a, Z[:, q, 3:]], -1)
    def init_fn(B):
        a = torch.randn(B, m, 3, device=dev); a = a / a.norm(dim=-1, keepdim=True)
        d = -1 + 0.3 * torch.randn(B, m, 1, device=dev)
        return torch.cat([a, d], -1)
    return init_fn, map_fn, reps

if __name__ == '__main__':
    import k12
    for S in k12.all59():
        gs = candidate_groups(S.name)
        if gs:
            print(S.name.split('_')[-1], [(lab, len(r)) for lab, r in gs])
