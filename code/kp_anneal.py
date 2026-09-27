"""Staged/annealed search: phase 1 explores broadly with an EQUAL-weight loss (like the original search.py) from
random starts, to find diverse basins with a low TOTAL violation count; phase 2 reseeds from a diverse pool of
those basins and switches to the heavy-X/light-S weighting to try to drive piercings to 0 while keeping few
self-crossings. Tries to escape the single strong (0,4,4) attractor found by seeding directly from the known
S4 near-miss."""
import time, sys
import numpy as np
import torch
import k12
import kp_search as kps


def anneal(S, t1=90.0, t2=90.0, B=4096, log=print, wX2=40.0, wS2=1.0, wConc2=0.2, noise=0.0):
    log(f'--- phase 1 (equal weight, random init) map {S.name} ---')
    top1, info1 = kps.run(S, B=B, time_limit=t1, wX=1.0, wS=1.0, lr=0.02, restart_every=300, log=log,
                           keep_top=30, noise=noise)
    pool = info1['pool']
    log(f'phase1 done: {info1["secs"]:.0f}s, starts {info1["starts"]}, pool best {pool[0][0] if pool else None}, '
        f'top best {top1[0][0] if top1 else None}')
    seeds = np.stack([r[1] for r in (top1 + pool)])
    log(f'--- phase 2 (heavy-X, seeded from {len(seeds)} diverse phase-1 basins) ---')
    init_fn = kps.seeded_init(seeds, noise_range=(0.01, 0.3))
    top2, info2 = kps.run(S, B=B, time_limit=t2, wX=wX2, wS=wS2, wConc=wConc2, lr=0.012, restart_every=300,
                           init_fn=init_fn, log=log, keep_top=15, noise=noise)
    log(f'phase2 done: {info2["secs"]:.0f}s, starts {info2["starts"]}, best {top2[0][0] if top2 else None}')
    return top1, top2, info1, info2


if __name__ == '__main__':
    num = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    t1 = float(sys.argv[2]) if len(sys.argv) > 2 else 90.0
    t2 = float(sys.argv[3]) if len(sys.argv) > 3 else 90.0
    S = k12.all59()[num - 1]
    top1, top2, info1, info2 = anneal(S, t1=t1, t2=t2)
    print('FINAL top2:')
    for k, _ in top2[:10]:
        print(k)
    np.save(f'kp_tmp_anneal_{num}.npy', np.stack([t[1] for t in top2]))
