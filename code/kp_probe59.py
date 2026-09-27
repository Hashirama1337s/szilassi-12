"""Broad probe: for each of the 59 maps, a short heavy-X/light-S weighted GPU search (random init, or seeded
from an existing near_*.npy snapshot when one exists for that map), to rank maps by best (piercings, faces
crossed, points) reached in a short budget. Used to decide where to spend the deeper search time."""
import glob, os, re, time, sys
import numpy as np
import torch
import k12
import kp_search as kps
import kp_common as kc

HERE = os.path.dirname(os.path.abspath(__file__))


def find_seed(num):
    pats = [f'near_A_{num}_*.npy', f'near_B_{num}_*.npy', f'near_C_{num}.npy', f'near_D_{num}.npy']
    for pat in pats:
        fs = glob.glob(os.path.join(HERE, pat))
        if fs:
            return fs[0]
    return None


def main(time_per_map=15.0, maps=None, out_log='kp_probe59_results.txt'):
    L = k12.all59()
    maps = maps or list(range(1, 60))
    results = []
    t00 = time.time()
    for num in maps:
        S = L[num - 1]
        seedf = find_seed(num)
        if seedf:
            base = np.load(seedf)
            init_fn = kps.seeded_init(base, noise_range=(0.03, 0.5))
            tag = os.path.basename(seedf)
        else:
            init_fn = kps.random_init(S, 'cuda')
            tag = 'random'
        try:
            top, info = kps.run(S, B=4096, time_limit=time_per_map, wX=40, wS=1.0, wConc=0.3, lr=0.015,
                                 restart_every=300, init_fn=init_fn)
            best = top[0][0] if top else (999, 999, 999)
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            best = (998, 998, 998)
        results.append((num, best, tag))
        print(f'map {num:2d}  seed={tag:20s}  best(pierce,faces,pts)={best}   t={time.time()-t00:6.0f}s', flush=True)
    results.sort(key=lambda x: x[1])
    with open(os.path.join(HERE, out_log), 'w') as f:
        for num, best, tag in results:
            f.write(f'{num}\t{best}\t{tag}\n')
    print('=== top 15 ===')
    for num, best, tag in results[:15]:
        print(num, best, tag)
    return results


if __name__ == '__main__':
    tpm = float(sys.argv[1]) if len(sys.argv) > 1 else 15.0
    main(time_per_map=tpm)
