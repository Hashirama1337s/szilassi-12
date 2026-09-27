"""Longer-patience unrestricted GPU runs on selected surfaces. usage: python gpu_long.py SECS k1,k2,..."""
import sys, json, time, numpy as np
from fractions import Fraction
import k12, search, check
secs, keys = float(sys.argv[1]), sys.argv[2].split(',')
surfs = {s.name.split('#')[-1]: s for s in k12.all59()}
out = open('sweep_D.log', 'a', encoding='utf-8'); out.write(f'=== {time.ctime()} long {secs}s {keys}\n'); out.flush()
for key in keys:
    S = surfs[key]
    found, info = search.run(S, B=4096, time_limit=secs, noise=0.005, restart_every=3000, lr=0.01, log=lambda s: None, seed=100 + int(key))
    res = []
    for H in found[:5]:
        Hq = [[Fraction(round(x * 10**7), 10**7) for x in h] for h in H]
        try: ok = check.check(S, Hq)[0]
        except ValueError: ok = 'degenerate'
        res.append(ok)
        if ok is True:
            json.dump(dict(surface=S.name, planes=[[str(x) for x in h] for h in Hq]), open(f'FOUND_{key}_long_{len(res)}.json', 'w'))
    np.save(f'near_D_{key}.npy', np.array([x[2] for x in info['near']]) if info['near'] else np.zeros(0))
    out.write(f'#{key:3s} long starts {info["starts"]:7d} min hard {info["minhard"]:3d} found {len(found)} exact {res} near {[x[0] for x in info["near"][:5]]}\n'); out.flush()
out.write('DONE\n')
