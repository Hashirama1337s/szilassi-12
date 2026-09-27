"""Search WITHIN a symmetric parametrisation (via symgroup.candidate_groups + kp_search's map_fn support): the
heavy-X/light-S weighted loss on the reduced d.o.f., aiming for 0 piercings with as few violated face-ORBITS as
possible (1 orbit = the fewest faces a symmetric config can ever cross, since violations come in whole orbits).
This is a natural way to beat the S4 near-miss's forced 4-face orbit: a C2 (order 2) or C3 (order 3) map only
forces 2 or 3 faces if a similarly-good 0-piercing point exists in its (much smaller) symmetric subspace."""
import time, sys
import numpy as np
import k12
import symgroup
import kp_search as kps


def sweep(labels=('C2', 'C3'), time_per=15.0, out='kp_sym_sweep_results.txt'):
    L = k12.all59()
    results = []
    t00 = time.time()
    for S in L:
        num = int(S.name.rsplit('_#', 1)[1])
        for lab, rho in symgroup.candidate_groups(S.name, n=S.n):
            if lab not in labels:
                continue
            init_fn, map_fn, reps = symgroup.param(rho, n=S.n, dev='cuda')
            top, info = kps.run(S, B=4096, time_limit=time_per, wX=40, wS=1.0, lr=0.015, restart_every=300,
                                 init_fn=init_fn, map_fn=map_fn, keep_top=5)
            best = top[0][0] if top else (999, 999, 999)
            results.append((num, lab, best))
            print(f'map {num:2d} {lab:3s}  best(pierce,faces,pts)={best}   t={time.time()-t00:6.0f}s', flush=True)
    results.sort(key=lambda x: x[2])
    with open(out, 'w') as f:
        for num, lab, best in results:
            f.write(f'{num}\t{lab}\t{best}\n')
    print('=== top 15 ===')
    for r in results[:15]:
        print(r)
    return results


if __name__ == '__main__':
    tp = float(sys.argv[1]) if len(sys.argv) > 1 else 15.0
    sweep(time_per=tp)
