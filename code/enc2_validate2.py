"""enc2_validate2 -- second validation batch: (a) first-violation configurations on rays leaving embedded Szilassi
realisations (7-vertex torus), (b) longer 12-face hill climbs (lower violation counts), (c) near-degenerate
configurations placed exactly at a violated condition of a near-miss (a plane moved to within 1e-30 of a point used
by that condition)."""
import sys, os, json, time
from fractions import Fraction as Fr
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k12
from enc2_core import surface
from enc2_geom import Geo
from enc2_validate import (validate, enc_for, quantize, rand_planes, climb, score, summary, GEO_FAMS)


def main():
    rng = np.random.default_rng(777)
    T7 = surface('t7')
    E7 = enc_for(T7)
    results = []
    t0 = time.time()

    def run(S, Hq, label):
        r = validate(S, Hq, label)
        results.append(r)
        if r['degenerate']:
            print(label, 'DEGENERATE', flush=True)
            return r
        print('%-34s %-4s emb=%d enc_ok=%d  viol(S/R/P/L)=%s  mism=%d  %.0fs' % (
            label, r['surface'], r['embedded'], r['enc_ok'],
            '/'.join(str(r['fam'][f][1]) for f in GEO_FAMS), r['n_mismatch'], time.time() - t0), flush=True)
        return r

    # (a) embedded torus realisations, then the first violation along random rays
    emb = []
    while len(emb) < 4:
        H, v = climb(E7, rng, rand_planes(rng, 7, 'sphere'), steps=600)
        if v == 0:
            emb.append(H)
    for k, H in enumerate(emb):
        for m in range(5):
            D = rng.normal(size=H.shape)
            lo, hi = 0.0, 1e-3
            while score(E7, (H + hi * D)[None])[0] == 0:
                lo, hi = hi, hi * 2
            for _ in range(30):                          # bisect to the first violation
                mid = (lo + hi) / 2
                if score(E7, (H + mid * D)[None])[0] == 0:
                    lo = mid
                else:
                    hi = mid
            run(T7, quantize(H + hi * D, 10 ** 9), 't7-firstviol-%d-%d' % (k, m))
            run(T7, quantize(H + lo * D, 10 ** 9), 't7-lastgood-%d-%d' % (k, m))
    # (b) longer 12-face climbs
    L = k12.all59()
    near = []
    for k in range(4):
        S = L[int(rng.integers(59))]
        E = enc_for(S)
        H, v = climb(E, rng, rand_planes(rng, 12, 'sphere'), steps=1500)
        Hq = quantize(H)
        r = run(S, Hq, 'climb12-long-%d(v=%d)' % (k, v))
        near.append((S, Hq, r))
    # (c) near-degenerate at a violated condition: move a plane of that condition onto the condition's point
    for S, Hq, r in near:
        E = enc_for(S)
        G = Geo(S, Hq)
        bad = [(f, key) for f in ('SIMP', 'PIERCE', 'REFL') for (ff, key, _) in E.conds if ff == f and G.cond_bad(f, key)]
        for t in range(3):
            fam, key = bad[int(rng.integers(len(bad)))]
            if fam == 'SIMP':
                i, s, t_ = key
                pt, l = (i, s, t_), E.nbr[i][s][0]           # vertex ist, moved onto plane (neighbour of s)
            elif fam == 'PIERCE':
                i, c, j = key
                pt, l = (i, j, c), E.nbr[j][c][0]
            else:
                tri, i, j = key
                pt = tri
                l = [x for x in E.S.link[i] if x not in tri][0]
            if l in pt:
                continue
            X = G.pt(*pt)
            a = Hq[l][:3]
            H2 = list(Hq)
            H2[l] = tuple(a) + (-sum(a[q] * X[q] for q in range(3)) + Fr(int(rng.choice([-1, 1])), 10 ** 30),)
            run(S, H2, 'at-boundary-%s-%s' % (fam, S.name.split('_')[-1]))
    return results


if __name__ == '__main__':
    import psutil
    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    res = main()
    tot = summary(res)
    with open(os.path.join(HERE, 'enc2_validate2_out.json'), 'w') as fh:
        json.dump(dict(results=res, totals=tot), fh, default=str, indent=1)
