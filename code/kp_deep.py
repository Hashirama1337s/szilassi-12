"""Deep-dive on one map: run kp_anneal's staged search, then for the best score reached with 0 (or few) genuine
piercings, try the keep-face reweighting refinement (bias the S loss to concentrate crossings onto a subset of
faces) starting from that result, over each candidate 'keep' face-set. Prints progress; does not write kp_best_
files itself (use kp_finalize.py on whatever .npy it saves once you've picked a winner)."""
import sys, itertools
import numpy as np
import torch
import k12
import kp_search as kps
import kp_anneal as kpa
import kp_common as kc


def refine_keep_faces(S, base_H, violated_faces, wX=40.0, time_limit=40.0, B=4096, heavy=8.0):
    results = {}
    for keep_n in range(0, min(3, len(violated_faces)) + 1):
        for keep in itertools.combinations(violated_faces, keep_n):
            fw = np.ones(S.n, dtype=np.float32)
            for f in violated_faces:
                if f not in keep:
                    fw[f] = heavy
            init_fn = kps.seeded_init(base_H, noise_range=(0.02, 0.3))
            top, info = kps.run(S, B=B, time_limit=time_limit, wX=wX, wS=1.0, face_weight=fw, init_fn=init_fn,
                                 lr=0.012, restart_every=300, keep_top=5)
            results[keep] = top
            print('  keep', keep, '->', top[0][0] if top else None)
    return results


def main(num, t1=90.0, t2=90.0, refine_time=40.0):
    S = k12.all59()[num - 1]
    top1, top2, info1, info2 = kpa.anneal(S, t1=t1, t2=t2)
    print('anneal best:', top2[0][0])
    best_key, best_H = top2[0]
    np.save(f'kp_tmp_deep_{num}_anneal.npy', best_H)
    if best_key[0] == 0 and best_key[1] > 0:
        # exact-check to get the true violated-face list (proxy key may differ slightly from exact checker)
        S2 = kc.surface_by_num(num)
        H = kc.to_rationals(best_H)
        ok, rep, gX, sc = kc.full_report(S2, H)
        print('exact score of anneal best:', sc)
        violated_faces = sorted({i for (i, c, d) in rep['S']})
        print('violated faces:', violated_faces)
        if sc[0] == 0 and len(violated_faces) > 1:
            results = refine_keep_faces(S, best_H, violated_faces, time_limit=refine_time)
            allrows = [r for rs in results.values() for r in rs]
            allrows.sort(key=lambda x: x[0])
            print('best after keep-face refine:', allrows[0][0])
            np.save(f'kp_tmp_deep_{num}_refined.npy', np.stack([r[1] for r in allrows[:10]]))
    return top2


if __name__ == '__main__':
    num = int(sys.argv[1])
    t1 = float(sys.argv[2]) if len(sys.argv) > 2 else 90.0
    t2 = float(sys.argv[3]) if len(sys.argv) > 3 else 90.0
    rt = float(sys.argv[4]) if len(sys.argv) > 4 else 40.0
    main(num, t1, t2, rt)
