"""Render a kp_best_<map>.json (or a raw npy row) as a matplotlib 3-D figure: 12 faces in distinct colours,
a normal view and an exploded view, with self-crossing points marked. Read-only w.r.t. existing modules."""
import sys, os, json, argparse, itertools
import numpy as np
from fractions import Fraction as Q
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import kp_common as kc
from check import point as _point, det3


def all_points_float(S, H):
    n = S.n
    P = {}
    for t in itertools.combinations(range(n), 3):
        if det3(H[t[0]][:3], H[t[1]][:3], H[t[2]][:3]) == 0:
            v = None
        else:
            v = _point(H, *t)
        for pm in itertools.permutations(t):
            P[pm] = v
    Pf = {k: (None if v is None else tuple(float(x) for x in v)) for k, v in P.items()}
    return Pf


def face_polygons(S, Pf):
    n = S.n
    m = n - 1
    faces = {}
    for i in range(n):
        L = S.link[i]
        faces[i] = [Pf[(i, L[s], L[(s + 1) % m])] for s in range(m)]
    return faces


def crossing_points(S, Pf, Sentries):
    pts = []
    for (i, c, d) in Sentries:
        pts.append((i, Pf[(i, c, d)]))
    return pts


COLORS = plt.cm.tab20(np.linspace(0, 1, 12))


def plot_one(ax, faces, cross_pts, explode=0.0, alpha=0.55):
    all_pts = np.array([p for poly in faces.values() for p in poly])
    centroid = all_pts.mean(0)
    for i, poly in faces.items():
        arr = np.array(poly)
        if explode:
            fc = arr.mean(0)
            direction = fc - centroid
            nrm = np.linalg.norm(direction)
            if nrm > 1e-9:
                arr = arr + explode * direction
        coll = Poly3DCollection([arr], alpha=alpha, facecolor=COLORS[i % 12], edgecolor='k', linewidths=0.6)
        ax.add_collection3d(coll)
    for (i, p) in cross_pts:
        pp = np.array(p)
        if explode:
            arr0 = np.array(faces[i])
            fc = arr0.mean(0)
            direction = fc - centroid
            nrm = np.linalg.norm(direction)
            if nrm > 1e-9:
                pp = pp + explode * direction
        ax.scatter(*pp, color='black', s=45, marker='o', depthshade=False, zorder=10)
    ax.set_box_aspect([1, 1, 1])
    lim = np.abs(all_pts).max() * (1 + explode * 1.5 + 0.1)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_axis_off()
    ax.view_init(elev=18, azim=35)


def render(mapnum, H, title, out_png, crossing_faces_pts):
    S = kc.surface_by_num(mapnum)
    Pf = all_points_float(S, H)
    faces = face_polygons(S, Pf)
    fig = plt.figure(figsize=(13, 6.5))
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    plot_one(ax1, faces, crossing_faces_pts, explode=0.0)
    ax1.set_title(f'{title}\n(assembled)')
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    plot_one(ax2, faces, crossing_faces_pts, explode=0.35)
    ax2.set_title('exploded (faces pulled outward)')
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    print('saved', out_png)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('map_num', type=int)
    ap.add_argument('path', help='.npy (row) or kp_best_*.json')
    ap.add_argument('--row', type=int, default=0)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    if args.path.endswith('.json'):
        data = json.load(open(args.path))
        H = [tuple(Q(x) for x in row) for row in data['planes']]
        sc = data['score']
        title = f"map #{args.map_num}: {sc['genuine_piercings']} piercings, {sc['faces_crossed']} faces crossed, {sc['crossing_points']} points"
        Sentries = [tuple(e) for e in data['report']['S']]
    else:
        arr = np.load(args.path)
        Hf = arr[args.row] if arr.ndim == 3 else arr
        H = kc.to_rationals(Hf)
        S = kc.surface_by_num(args.map_num)
        ok, rep, gX, sc = kc.full_report(S, H)
        title = f"map #{args.map_num}: {sc[0]} piercings, {sc[1]} faces crossed, {sc[2]} points"
        Sentries = rep['S']
    S = kc.surface_by_num(args.map_num)
    Pf = all_points_float(S, H)
    cpts = crossing_points(S, Pf, Sentries)
    out_png = args.out or f'kp_best_{args.map_num}.png'
    render(args.map_num, H, title, out_png, cpts)


if __name__ == '__main__':
    main()
