"""Finalize a candidate: exact-verify (checker A + genuine-piercing filter), degeneracy check, and if it's a
genuine improvement, save kp_best_<map>.json + kp_best_<map>.png and append a line to kp_log.txt.
Usage: python kp_finalize.py <map_num> <path.npy> [--row R] [--label L] [--force]"""
import argparse, os, json, subprocess, sys
import numpy as np
import kp_common as kc

HERE = os.path.dirname(os.path.abspath(__file__))


def current_best(mapnum):
    fn = os.path.join(HERE, f'kp_best_{mapnum}.json')
    if not os.path.exists(fn):
        return None
    d = json.load(open(fn))
    s = d['score']
    return (s['genuine_piercings'], s['faces_crossed'], s['crossing_points'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('map_num', type=int)
    ap.add_argument('path')
    ap.add_argument('--row', type=int, default=0)
    ap.add_argument('--label', default='')
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--grid', type=int, default=10_000_000)
    args = ap.parse_args()

    arr = np.load(args.path)
    Hf = arr[args.row] if arr.ndim == 3 else arr
    S = kc.surface_by_num(args.map_num)
    H = kc.to_rationals(Hf, args.grid)
    ok, rep, gX, sc = kc.full_report(S, H)
    mer, rms = kc.min_edge_ratio(S, H)
    print('map', args.map_num, 'score', sc, 'embedded', ok, 'min_edge_ratio', f'{mer:.5f}')

    cb = current_best(args.map_num)
    is_improve = args.force or cb is None or sc < cb
    if not is_improve:
        print(f'not an improvement over existing kp_best_{args.map_num}.json (score {cb}); not saved (use --force).')
        return

    degen_note = ''
    if mer < 0.01:
        degen_note = f' DEGENERATE WARNING: min_edge/rms={mer:.5f} < 0.01'

    out = {
        'map': f'manifold_lex_d2_n12_o1_g6_#{args.map_num}',
        'score': {'genuine_piercings': sc[0], 'faces_crossed': sc[1], 'crossing_points': sc[2]},
        'embedded': ok,
        'min_edge_over_rms_radius': mer,
        'planes': [[str(x) for x in h] for h in H],
        'report': {'S': rep['S'], 'X_raw': rep['X'], 'X_genuine': gX, 'G': rep['G']},
        'source': {'path': args.path, 'row': args.row, 'label': args.label},
    }
    jf = os.path.join(HERE, f'kp_best_{args.map_num}.json')
    with open(jf, 'w') as f:
        json.dump(out, f, indent=1)
    print('saved', jf)

    pf = os.path.join(HERE, f'kp_best_{args.map_num}.png')
    subprocess.run([sys.executable, os.path.join(HERE, 'kp_render.py'), str(args.map_num), jf, '--out', pf],
                    check=True)

    kc.log_line(f'map #{args.map_num}  score={sc}  embedded={ok}  min_edge_ratio={mer:.5f}  '
                f'label={args.label}{degen_note}')


if __name__ == '__main__':
    main()
