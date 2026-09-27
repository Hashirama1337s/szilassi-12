"""Standalone exact verification of a candidate configuration: rounds to a rational grid, runs checker A
(check.py) plus the genuine-piercing filter, reports the lexicographic score, degeneracy (min edge length /
RMS radius, min |det| at a vertex), and optionally saves a kp_best_<map>.json.
Usage: python kp_verify.py <map_num> <path.npy> [--row R] [--save]"""
import sys, os, json, argparse
import numpy as np
import kp_common as kc


def verify(num, H_float, grid=10_000_000):
    S = kc.surface_by_num(num)
    H = kc.to_rationals(H_float, grid)
    ok, rep, gX, sc = kc.full_report(S, H)
    mer, rms = kc.min_edge_ratio(S, H)
    return S, H, ok, rep, gX, sc, mer, rms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('map_num', type=int)
    ap.add_argument('path')
    ap.add_argument('--row', type=int, default=0)
    ap.add_argument('--grid', type=int, default=10_000_000)
    ap.add_argument('--save', action='store_true')
    args = ap.parse_args()

    arr = np.load(args.path)
    if arr.ndim == 3:
        Hf = arr[args.row]
    else:
        Hf = arr
    S, H, ok, rep, gX, sc, mer, rms = verify(args.map_num, Hf, args.grid)
    print('map', args.map_num, 'embedded' if ok else 'NOT embedded (near-miss)')
    print('score (genuine_pierce, faces_crossed, crossing_points) =', sc)
    print('raw S entries (i,c,d):', rep['S'])
    print('raw X entries (i,j,k):', rep['X'], '  genuine X:', gX)
    print('G (genericity issues):', rep['G'])
    print(f'min edge / rms_radius = {mer:.5f}   rms_radius = {rms:.4f}')
    if args.save:
        out = {
            'map': f'manifold_lex_d2_n12_o1_g6_#{args.map_num}',
            'score': {'genuine_piercings': sc[0], 'faces_crossed': sc[1], 'crossing_points': sc[2]},
            'embedded': ok,
            'min_edge_over_rms': mer,
            'planes': [[str(x) for x in h] for h in H],
            'report': {'S': rep['S'], 'X_raw': rep['X'], 'X_genuine': gX, 'G': rep['G']},
        }
        fn = f'kp_best_{args.map_num}.json'
        with open(fn, 'w') as f:
            json.dump(out, f, indent=1)
        print('saved', fn)


if __name__ == '__main__':
    main()
