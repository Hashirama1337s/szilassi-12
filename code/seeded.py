"""Seeded unrestricted search: start from near-miss configurations (e.g. symmetric ones), perturb, optimise all 24 dof."""
import sys, json, numpy as np, torch
from fractions import Fraction
import k12, search, check
key, npy, secs, sig = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
S = {s.name.split('_')[-1]: s for s in k12.all59()}['#' + key]
seeds = torch.tensor(np.load(npy)[:5], dtype=torch.float32, device='cuda')
def init_fn(B):
    q = torch.randint(0, len(seeds), (B,), device='cuda')
    return seeds[q] + sig * torch.randn(B, S.n, 4, device='cuda')
found, info = search.run(S, B=2048, time_limit=secs, noise=0.002, param=(init_fn, lambda Z: Z), log=lambda s: None, seed=1,
                         margin=0.005)
res = []
for H in found[:5]:
    Hq = [[Fraction(round(x * 10**7), 10**7) for x in h] for h in H]
    try: ok = check.check(S, Hq)[0]
    except ValueError as e: ok = str(e)
    res.append(ok)
    if ok is True:
        json.dump(dict(surface=S.name, planes=[[str(x) for x in h] for h in Hq]), open(f'FOUND_{key}_seeded_{len(res)}.json', 'w'))
print(f'#{key} seeded sig {sig}: starts {info["starts"]} min hard {info["minhard"]} found {len(found)} exact {res} near {[x[0] for x in info["near"][:6]]}', flush=True)
np.save(f'near_seeded_{key}.npy', np.array([x[2] for x in info['near']]) if info['near'] else np.zeros(0))
