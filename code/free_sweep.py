"""Unrestricted numerical sweep over all 59 surfaces. usage: python free_sweep.py SECS TAG [first last]"""
import sys, json, time, numpy as np
from fractions import Fraction
import k12, search, check
secs, tag = float(sys.argv[1]), sys.argv[2]
lo, hi = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (1, 59)
out = open(f'sweep_{tag}.log', 'a', encoding='utf-8'); out.write(f'=== {time.ctime()} free {secs}s #{lo}-#{hi}\n'); out.flush()
for S in k12.all59()[lo - 1:hi]:
    key = S.name.split('#')[-1]
    found, info = search.run(S, B=4096, time_limit=secs, noise=0.003, log=lambda s: None, seed=int(key))
    res = []
    for H in found[:5]:
        Hq = [[Fraction(round(x * 10**7), 10**7) for x in h] for h in H]
        try: ok = check.check(S, Hq)[0]
        except ValueError as e: ok = 'degenerate'
        res.append(ok)
        if ok is True:
            json.dump(dict(surface=S.name, planes=[[str(x) for x in h] for h in Hq]), open(f'FOUND_{key}_free_{len(res)}.json', 'w'))
    np.save(f'near_{tag}_{key}.npy', np.array([x[2] for x in info['near']]) if info['near'] else np.zeros(0))
    out.write(f'#{key:3s} free starts {info["starts"]:8d} min hard {info["minhard"]:3d} found {len(found)} exact {res} near {[x[0] for x in info["near"][:5]]}\n'); out.flush()
out.write('DONE\n')
