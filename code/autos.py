"""Automorphism groups of the 59 surfaces (and the 7-torus): each automorphism is fixed by the image of one flag."""
import k12, itertools, json, collections

def automorphisms(S):
    tri = set(S.tris)
    otri = set()
    for (a, b, c) in S.otris:
        otri |= {(a, b, c), (b, c, a), (c, a, b)}
    a0, b0, c0 = S.otris[0]
    out = []
    for t in S.tris:
        for img in itertools.permutations(t):
            f = {a0: img[0], b0: img[1], c0: img[2]}
            todo = [(a0, b0, c0)]
            ok = True
            seen = {tuple(sorted((a0, b0, c0)))}
            while todo and ok:
                x, y, z = todo.pop()
                for (p, q, r) in ((x, y, z), (y, z, x), (z, x, y)):
                    e = tuple(sorted((p, q)))
                    w = [u for u in S.elink[e] if u != r][0]
                    fw = [u for u in S.elink[tuple(sorted((f[p], f[q])))] if u != f[r]][0]
                    if w in f:
                        if f[w] != fw: ok = False; break
                    else:
                        if fw in f.values(): ok = False; break
                        f[w] = fw
                    s = tuple(sorted((p, q, w)))
                    if s not in seen:
                        seen.add(s); todo.append((q, p, w))
            if ok and len(f) == S.n and all(tuple(sorted(f[v] for v in t2)) in tri for t2 in S.tris):
                perm = tuple(f[v] for v in range(S.n))
                pres = (perm[a0], perm[b0], perm[c0]) in otri
                out.append((perm, pres))
    return out

def order(p):
    q, k = p, 1
    while q != tuple(range(len(p))):
        q = tuple(p[x] for x in q); k += 1
    return k

res = {}
for idx, S in enumerate([k12.Surface(k12.torus7(), 't7')] + k12.all59()):
    A = automorphisms(S)
    stats = collections.Counter((order(p), pres, sum(p[v] == v for v in range(S.n))) for p, pres in A)
    res[S.name] = dict(size=len(A), pres=sum(1 for _, pr in A if pr), elems=[[list(p), pr, order(p)] for p, pr in A if order(p) > 1])
    print(idx, S.name.split('_')[-1], '|Aut|', len(A), '|Aut+|', res[S.name]['pres'], sorted((o, pr, fx) for (o, pr, fx) in stats if o > 1))
json.dump(res, open('automorphisms.json', 'w'))
