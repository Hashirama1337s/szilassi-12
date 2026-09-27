"""enc2_validate -- compare the enc2 clauses, evaluated on the exact chirotope of a configuration, with the exact
geometric truth (enc2_geom), condition by condition, plus the overall embedding test.

Configurations: random planes (two distributions), near-misses from a hill climb on the number of violated
conditions (free and symmetric), embedded Szilassi realisations found by the same climb on the 7-vertex torus and their
perturbations, and near-degenerate configurations (a plane moved to pass within 1e-30 of a dual vertex)."""
import sys, os, time, itertools, json, math
from fractions import Fraction as Fr
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k12
from enc2_core import Enc, surface, perm_sign
from enc2_geom import Geo

GEO_FAMS = ('SIMP', 'REFL', 'PIERCE', 'LOC')
CLIMB_FAMS = ('SIMP', 'REFL', 'PIERCE')

_ENC = {}


def enc_for(S):
    if S.name not in _ENC:
        _ENC[S.name] = Enc(S)
    return _ENC[S.name]


def quantize(H, scale=10 ** 6):
    return [tuple(Fr(int(round(float(x) * scale))) for x in h) for h in H]


def validate(S, Hq, label, sym=None):
    """Hq: exact planes.  Returns a result dict."""
    E = enc_for(S)
    G = Geo(S, Hq)
    chi = G.chi_exact(E.bases)
    if any(c == 0 for c in chi):
        return dict(label=label, surface=S.name, degenerate=True)
    V = E.evaluate(np.array(chi, dtype=np.int64)[:, None] > 0)
    assert E.defs_ok(V), 'Tseitin definitions violated'
    res = dict(label=label, surface=S.name.split('_')[-1], degenerate=False, fam={})
    allfams = ('GP',) + GEO_FAMS
    idx, viol = E.violated(V, allfams)
    viol = viol[:, 0]
    mism = []
    for k, ci in enumerate(idx):
        fam, key, _ = E.conds[ci]
        enc_bad = bool(viol[k])
        f = res['fam'].setdefault(fam, [0, 0, 0])        # compared, violated(geometry), mismatches
        f[0] += 1
        if fam == 'GP':
            geo_bad = False                             # GP holds for every real chirotope
        else:
            geo_bad = G.cond_bad(fam, key)
        f[1] += geo_bad
        if geo_bad != enc_bad:
            f[2] += 1
            mism.append((fam, key, geo_bad, enc_bad))
    enc_all_ok = not viol.any()
    emb, nsimp, npair = G.embedded(detail=True)
    res.update(embedded=emb, enc_ok=enc_all_ok, nonsimple_faces=nsimp, bad_pairs=npair, mismatches=mism[:10],
               n_mismatch=len(mism), overall_match=(emb == enc_all_ok))
    if sym is not None:                                 # check the symmetric substitution on this exact config
        Es = sym
        ok = True
        for v in range(1, E.nchi + 1):
            r = Es.rep[v]
            if (chi[v - 1] > 0) != ((chi[abs(r) - 1] > 0) == (r > 0)):
                ok = False
        res['sym_subst_ok'] = ok
    return res


# ---------------- generators ----------------
def rand_planes(rng, n, kind):
    if kind == 'gauss':
        return np.concatenate([rng.normal(size=(n, 3)), rng.normal(size=(n, 1))], 1)
    a = rng.normal(size=(n, 3))
    a /= np.linalg.norm(a, axis=1, keepdims=True)
    return np.concatenate([a, -1 + 0.3 * rng.normal(size=(n, 1))], 1)


def score(E, Hb):
    chi, mind = E.chi_from_H(Hb)
    V = E.evaluate(chi)
    _, viol = E.violated(V, CLIMB_FAMS)
    s = viol.sum(axis=0).astype(float)
    s[mind < 1e-9] = 1e9
    return s


def climb(E, rng, H, steps, batch=48, target=0, param=None, sigma=0.25, verbose=False):
    """(1+lambda) hill climb on the number of violated SIMP/REFL/PIERCE conditions.
    param: None (free planes) or (reps, face_src) for symmetric parametrisation."""
    if param is None:
        mk = lambda Z: Z
    else:
        mk = param
    cur = score(E, mk(H)[None])[0]
    for st in range(steps):
        C = np.repeat(H[None], batch, 0)
        w = rng.integers(H.shape[0], size=batch)
        sig = sigma * np.exp(rng.uniform(-3, 0, size=batch))
        C[np.arange(batch), w] += sig[:, None] * rng.normal(size=(batch, 4))
        sc = score(E, np.stack([mk(c) for c in C]))
        b = int(np.argmin(sc))
        if sc[b] <= cur:
            H, cur = C[b], sc[b]
        if verbose and st % 50 == 0:
            print('   step', st, 'viol', cur, flush=True)
        if cur <= target:
            break
    return H, cur


def sym_param(rho, n):
    reps, src = [], {}
    for f in range(n):
        if f in src:
            continue
        reps.append(f)
        for g, M in rho.items():
            src[g[f]] = (len(reps) - 1, np.asarray(M, float))
    def mk(Z):
        out = np.zeros((n, 4))
        for f in range(n):
            r, M = src[f]
            out[f, :3] = M @ Z[r, :3]
            out[f, 3] = Z[r, 3]
        return out
    return reps, src, mk


def rational_rep(rho):
    """conjugate rho into SO(3)-equivalent rational matrices (entries 0,+-1): z-axis -> (1,1,1), x-axis -> (1,-1,0)."""
    if all(np.allclose(M, np.round(M), atol=1e-12) for M in rho.values()):
        return {g: np.round(M).astype(int) for g, M in rho.items()}
    ex = np.array([1, -1, 0]) / math.sqrt(2)
    ez = np.array([1, 1, 1]) / math.sqrt(3)
    ey = np.cross(ez, ex)
    R = np.stack([ex, ey, ez], 1)                     # columns: images of x, y, z
    assert abs(np.linalg.det(R) - 1) < 1e-12
    out = {}
    for g, M in rho.items():
        Mp = R @ M @ R.T
        assert np.allclose(Mp, np.round(Mp), atol=1e-9)
        out[g] = np.round(Mp).astype(int)
    return out


def exact_sym_planes(Zq, src_rat, n):
    """Zq: orbit-rep planes (Fractions); src_rat: f -> (rep index, integer matrix)."""
    out = []
    for f in range(n):
        r, M = src_rat[f]
        a = [sum(int(M[x][y]) * Zq[r][y] for y in range(3)) for x in range(3)]
        out.append(tuple(a) + (Zq[r][3],))
    return out


def near_degenerate(G, Hq, rng):
    """move one plane l so that it passes within 1e-30 of a dual vertex (a triangle ijk not containing l)."""
    S = G.S
    tri = S.tris[rng.integers(len(S.tris))]
    l = int(rng.choice([x for x in range(S.n) if x not in tri]))
    X = G.pt(*tri)
    a = Hq[l][:3]
    eps = Fr(int(rng.choice([-1, 1])), 10 ** 30)
    d = -sum(a[k] * X[k] for k in range(3)) + eps
    H2 = list(Hq)
    H2[l] = tuple(a) + (d,)
    return H2, (tri, l)


def main(quick=False):
    import symgroup
    rng = np.random.default_rng(20260927)
    L = k12.all59()
    names = {S.name.split('_')[-1]: S for S in L}
    T7 = surface('t7')
    results = []
    t0 = time.time()

    def run(S, Hq, label, sym=None):
        r = validate(S, Hq, label, sym=sym)
        results.append(r)
        if r['degenerate']:
            print(label, 'DEGENERATE (skipped)', flush=True)
            return r
        print('%-34s %-4s emb=%d enc_ok=%d  viol(S/R/P/L)=%s  mism=%d  %s%.0fs' % (
            label, r['surface'], r['embedded'], r['enc_ok'],
            '/'.join(str(r['fam'][f][1]) for f in GEO_FAMS), r['n_mismatch'],
            ('symOK=%d ' % r['sym_subst_ok']) if 'sym_subst_ok' in r else '', time.time() - t0), flush=True)
        return r

    # 1. random configurations (12 faces, random surfaces), two distributions
    nr = 3 if quick else 12
    for k in range(nr):
        S = L[int(rng.integers(59))]
        kind = 'gauss' if k % 2 == 0 else 'sphere'
        run(S, quantize(rand_planes(rng, 12, kind), 1000 if kind == 'gauss' else 10 ** 6), 'random-%s-%d' % (kind, k))
    if quick:
        return results
    # 2. torus7: climb to embedded Szilassi realisations, then perturb (near-misses), + random
    t7_emb = []
    E7 = enc_for(T7)
    for k in range(40):
        H, v = climb(E7, rng, rand_planes(rng, 7, 'sphere'), steps=600)
        if v == 0:
            t7_emb.append(H)
            if len(t7_emb) >= 4:
                break
    print('torus7 embedded found:', len(t7_emb), flush=True)
    for k, H in enumerate(t7_emb):
        run(T7, quantize(H), 't7-embedded-%d' % k)
        for m in range(2):
            for tries in range(200):                    # perturb until 1..6 violations
                H2 = H + 0.02 * (1 + tries / 20) * rng.normal(size=H.shape)
                v = score(E7, H2[None])[0]
                if 1 <= v <= 6:
                    break
            run(T7, quantize(H2), 't7-nearmiss-%d-%d(v=%d)' % (k, m, v))
    for k in range(2):
        run(T7, quantize(rand_planes(rng, 7, 'gauss'), 1000), 't7-random-%d' % k)
    # 3. 12-face near-misses by free hill climbing
    for k in range(8):
        S = L[int(rng.integers(59))]
        E = enc_for(S)
        H, v = climb(E, rng, rand_planes(rng, 12, 'sphere'), steps=250)
        run(S, quantize(H), 'climb12-%d(v=%d)' % (k, v))
    # 4. exact symmetric near-misses (rational representations), check symmetric substitution too
    for key, lab in [('#54', 'T'), ('#59', 'D2'), ('#5', 'S4'), ('#45', 'C3'), ('#3', 'D3'), ('#48', 'D3')]:
        S = names[key]
        rho = dict(symgroup.candidate_groups(S.name))[lab]
        rr = rational_rep(rho)
        for g in rr:
            assert round(np.linalg.det(rr[g])) == round(np.linalg.det(rho[g]))
        reps, src, mk = sym_param({g: rr[g].astype(float) for g in rr}, 12)
        E = enc_for(S)
        Es = Enc(S, sym=[(g, int(round(np.linalg.det(M)))) for g, M in rho.items()])
        Z0 = rand_planes(rng, len(reps), 'sphere')
        Z, v = climb(E, rng, Z0, steps=200, param=mk)
        Zq = quantize(Z)
        src_rat = {f: (r, rr_M) for f, (r, rr_M) in ((f, (src[f][0], np.round(src[f][1]).astype(int))) for f in src)}
        Hq = exact_sym_planes(Zq, src_rat, 12)
        run(S, Hq, 'sym-%s-%s(v=%d)' % (key, lab, v), sym=Es)
    # 5. near-degenerate configurations
    base = [r for r in results if not r['degenerate']]
    pool = [('t7', T7), ('#54', names['#54']), ('#3', names['#3'])]
    for k in range(8):
        key, S = pool[k % 3]
        n = S.n
        Hq = quantize(rand_planes(rng, n, 'sphere'))
        G = Geo(S, Hq)
        H2, info = near_degenerate(G, Hq, rng)
        run(S, H2, 'near-degenerate-%d %s' % (k, key))
    return results


def summary(results):
    tot = {}
    nd = [r for r in results if not r['degenerate']]
    for r in nd:
        for f, (c, v, m) in r['fam'].items():
            t = tot.setdefault(f, [0, 0, 0])
            t[0] += c; t[1] += v; t[2] += m
    print('configs', len(nd), '(degenerate skipped: %d)' % (len(results) - len(nd)))
    print('family: compared / violated / mismatches')
    for f, (c, v, m) in tot.items():
        print('  %-7s %8d %8d %4d' % (f, c, v, m))
    print('overall embedded-vs-clauses mismatches:', sum(not r['overall_match'] for r in nd),
          '; embedded configs:', sum(r['embedded'] for r in nd))
    if any('sym_subst_ok' in r for r in nd):
        print('symmetric substitution consistent on exact symmetric configs:',
              all(r['sym_subst_ok'] for r in nd if 'sym_subst_ok' in r))
    for r in nd:
        if r['n_mismatch']:
            print('MISMATCH', r['label'], r['mismatches'])
    return tot


if __name__ == '__main__':
    quick = len(sys.argv) > 1 and sys.argv[1] == 'quick'
    res = main(quick=quick)
    tot = summary(res)
    with open(os.path.join(HERE, 'enc2_validate_out.json'), 'w') as fh:
        json.dump(dict(results=res, totals=tot), fh, default=str, indent=1)
