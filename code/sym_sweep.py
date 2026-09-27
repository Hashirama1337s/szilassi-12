"""Symmetric sweep over the 59 surfaces, smallest search spaces first. Every hit is rounded to rationals and exact-checked."""
import sys, json, time, os, numpy as np, torch
from fractions import Fraction
import k12, search, symgroup, check

surfs = {S.name.split('_')[-1]: S for S in k12.all59()}
labels = sys.argv[1].split(',')              # e.g. T,D3,D2,S4,C3
secs = float(sys.argv[2])
tag = sys.argv[3]
out = open(f'sweep_{tag}.log', 'a', encoding='utf-8')
def log(s):
    out.write(s + '\n'); out.flush()
jobs = []
for key, S in surfs.items():
    seen_lab = set()
    for lab, rho in symgroup.candidate_groups(S.name):
        if lab not in labels: continue
        if lab in ('C2', 'C3') and key != '#59' and lab in seen_lab: continue   # conjugate copies
        seen_lab.add(lab)
        jobs.append((12 // len(rho) * 4, key, lab, rho))
jobs.sort(key=lambda j: j[0])
log(f'=== {time.ctime()} labels {labels} secs {secs} jobs {len(jobs)}')
for npar, key, lab, rho in jobs:
    S = surfs[key]
    init_fn, map_fn, reps = symgroup.param(rho)
    found, info = search.run(S, B=4096, time_limit=secs, noise=0.003, param=(init_fn, map_fn), log=lambda s: None, seed=hash((key, lab)) % 1000)
    res = []
    for H in found[:5]:
        Hq = [[Fraction(round(x * 10**7), 10**7) for x in h] for h in H]
        try:
            ok, rep = check.check(S, Hq)
        except ValueError as e:
            ok, rep = False, str(e)
        res.append(ok)
        if ok:
            json.dump(dict(surface=S.name, group=lab, planes=[[str(x) for x in h] for h in Hq]), open(f'FOUND_{key[1:]}_{lab}_{len(res)}.json', 'w'))
    np.save(f'near_{tag}_{key[1:]}_{lab}.npy', np.array([x[2] for x in info['near']]) if info['near'] else np.zeros(0))
    log(f'{key:4s} {lab:3s} params {npar:2d}  starts {info["starts"]:8d}  min hard {info["minhard"]:3d}  found {len(found)}  exact {res}'
        f'  near {[x[0] for x in info["near"][:5]]}')
log('DONE')
