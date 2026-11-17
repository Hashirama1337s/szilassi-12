"""Figures for the szilassi-12 repository.

Usage:  python figures/make_figures.py <repo root>
Writes figures/*.svg next to this script. Inputs, all shipped in the repository:
  convex7/core.py            T3V, T3F: Szilassi's coordinates (Grunbaum & Szilassi 2009, Table 3), read as data with
                             the ast module (the module itself is not imported, so nothing is run or written there)
  results/certificates.csv   Part B: one row per proved case (map, group, solve_s, verdict)
Planarity, face adjacency, reflex corners and the convex face are recomputed here in exact rational arithmetic.
The 3-D picture is an orthographic view with exact per-pixel visibility (one ray per pixel, nearest face wins),
embedded as an image; it is an illustration of the exact coordinates, not part of any proof.
Checked copy (figure referee, 2026-11-17): the "un-dented face" label moved clear of face 5 in the view from below.
"""
import ast
import csv
import itertools
import os
import sys
from collections import Counter
from fractions import Fraction

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.patches import FancyBboxPatch, Polygon

BG, BORDER, GRID, AXC, TXT, ACC = "#1c2026", "#2f353d", "#22272e", "#9aa4ae", "#e8edf2", "#39ff14"
CAT = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c", "#9085e9", "#e66767"]


def style():
    plt.rcParams.update({
        "svg.fonttype": "none", "svg.hashsalt": "szilassi-12-figures",
        "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14, "xtick.labelsize": 12.5, "ytick.labelsize": 12.5,
        "legend.fontsize": 12.5, "figure.facecolor": "none", "axes.facecolor": BG, "axes.edgecolor": AXC,
        "axes.labelcolor": TXT, "axes.titlecolor": TXT, "xtick.color": AXC, "ytick.color": AXC, "text.color": TXT,
        "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0,
    })


def frame(fig):
    w, h = fig.get_size_inches()
    fig.patches.append(FancyBboxPatch((0.03, 0.03), w - 0.06, h - 0.06, boxstyle="round,pad=0,rounding_size=0.18",
                                      transform=fig.dpi_scale_trans, fc=BG, ec=BORDER, lw=1.5, zorder=-100))


def save(fig, name):
    here = os.path.dirname(os.path.abspath(__file__))
    out = here if os.path.basename(here) == "figures" else os.path.join(here, "figures")   # repo: figures/make_figures.py
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, name), format="svg", dpi=200, metadata={"Date": None, "Creator": None})
    plt.close(fig)
    print("wrote figures/" + name)


# ------------------------------------------------------------------ exact data
def read_table3(root):
    src = open(os.path.join(root, "convex7", "core.py"), encoding="utf-8").read()
    env = {"Q": Fraction}
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and getattr(node.targets[0], "id", "") in ("T3V", "T3F"):
            for sub in ast.walk(node.value):       # only literals, tuples, dicts and Q('...') calls are allowed
                assert isinstance(sub, (ast.Dict, ast.Tuple, ast.List, ast.Constant, ast.UnaryOp, ast.USub, ast.Call,
                                        ast.Name, ast.Load)), type(sub)
                if isinstance(sub, ast.Call):
                    assert sub.func.id == "Q"
                if isinstance(sub, ast.Name):
                    assert sub.id == "Q"
            env[node.targets[0].id] = eval(compile(ast.Expression(node.value), "core.py", "eval"), {"__builtins__": {}}, env)
    V = {k: tuple(Fraction(x) for x in v) for k, v in env["T3V"].items()}
    return V, [list(f) for f in env["T3F"]]


def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def cross(a, b): return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
def dot(a, b): return sum(x * y for x, y in zip(a, b))


def analyse(V, F):
    """Exact: plane normal per face, planarity, reflex corners, convex faces, shared edges."""
    normals, reflex = [], []
    for f in F:
        P = [V[v] for v in f]
        n = None
        for a, b in itertools.combinations(range(1, 6), 2):
            n = cross(sub(P[a], P[0]), sub(P[b], P[0]))
            if any(n):
                break
        assert all(dot(n, sub(p, P[0])) == 0 for p in P), "face not planar"
        area = sum(dot(cross(P[s], P[(s + 1) % 6]), n) for s in range(6))
        r = []
        for s in range(6):
            t = dot(cross(sub(P[s], P[s - 1]), sub(P[(s + 1) % 6], P[s])), n)
            assert t != 0
            r.append((t > 0) != (area > 0))
        normals.append(n)
        reflex.append(r)
    E = {}
    for i, f in enumerate(F):
        for s in range(6):
            E.setdefault(frozenset((f[s], f[(s + 1) % 6])), []).append(i)
    assert len(E) == 21 and all(len(x) == 2 for x in E.values())
    pairs = Counter(frozenset(x) for x in E.values())
    assert len(pairs) == 21 and set(pairs.values()) == {1}, "not every pair of faces shares exactly one edge"
    return normals, reflex, E


# ------------------------------------------------------------------ drawing helpers
def ear_clip(P2):
    idx = list(range(len(P2)))
    area = sum(P2[i][0] * P2[(i + 1) % len(P2)][1] - P2[(i + 1) % len(P2)][0] * P2[i][1] for i in idx)
    sgn = 1 if area > 0 else -1
    tris = []
    while len(idx) > 3:
        for k in range(len(idx)):
            a, b, c = idx[k - 1], idx[k], idx[(k + 1) % len(idx)]
            A, B, C = P2[a], P2[b], P2[c]
            if sgn * ((B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0])) <= 0:
                continue
            inside = False
            for j in idx:
                if j in (a, b, c):
                    continue
                Pj = P2[j]
                d1 = (B[0] - A[0]) * (Pj[1] - A[1]) - (B[1] - A[1]) * (Pj[0] - A[0])
                d2 = (C[0] - B[0]) * (Pj[1] - B[1]) - (C[1] - B[1]) * (Pj[0] - B[0])
                d3 = (A[0] - C[0]) * (Pj[1] - C[1]) - (A[1] - C[1]) * (Pj[0] - C[0])
                if sgn * d1 >= 0 and sgn * d2 >= 0 and sgn * d3 >= 0:
                    inside = True
                    break
            if not inside:
                tris.append((a, b, c))
                idx.pop(k)
                break
        else:
            raise RuntimeError("ear clipping failed")
    tris.append(tuple(idx))
    return tris


def plane_basis(n):
    n = np.array([float(x) for x in n])
    n /= np.linalg.norm(n)
    e1 = np.cross(n, [1.0, 0, 0]) if abs(n[0]) < 0.9 else np.cross(n, [0, 1.0, 0])
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(n, e1), n


def rotation(az, el):
    az, el = np.radians(az), np.radians(el)
    Rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(el), -np.sin(el)], [0, np.sin(el), np.cos(el)]])
    return Rx @ Rz          # view along +y after rotation: screen = (x, z), depth = y


def face_triangles(V, F, normals):
    out = []
    for i, f in enumerate(F):
        e1, e2, _ = plane_basis(normals[i])
        P3 = [np.array([float(x) for x in V[v]]) for v in f]
        P2 = [(p @ e1, p @ e2) for p in P3]
        T = ear_clip(P2)
        area = abs(sum(P2[k][0] * P2[(k + 1) % 6][1] - P2[(k + 1) % 6][0] * P2[k][1] for k in range(6))) / 2
        tri_area = sum(abs((P2[b][0] - P2[a][0]) * (P2[c][1] - P2[a][1]) - (P2[c][0] - P2[a][0]) * (P2[b][1] - P2[a][1]))
                       for a, b, c in T) / 2
        assert len(T) == 4 and abs(area - tri_area) < 1e-9 * area, "triangulation does not tile the face"
        for (a, b, c) in T:
            out.append((i, P3[a], P3[b], P3[c]))
    return out


def raycast(V, F, normals, az, el, npx=900, pad=2.0):
    """Exact-visibility image: for every pixel, the face hit first by a parallel ray along the view direction
    (orthographic). Returns the face index per pixel (-1 = background) and the screen extent."""
    R = rotation(az, el)
    tri = face_triangles(V, F, normals)
    allp = np.array([R @ np.array([float(x) for x in v]) for v in V.values()])
    x0, x1 = allp[:, 0].min() - pad, allp[:, 0].max() + pad
    z0, z1 = allp[:, 2].min() - pad, allp[:, 2].max() + pad
    s = max(x1 - x0, z1 - z0) / npx
    xs = np.arange(x0, x1, s) + s / 2
    zs = np.arange(z0, z1, s) + s / 2
    X, Z = np.meshgrid(xs, zs)
    best = np.full(X.shape, np.inf)
    fid = np.full(X.shape, -1)
    for (i, A, B, C) in tri:
        A, B, C = R @ A, R @ B, R @ C
        d = (B[0] - A[0]) * (C[2] - A[2]) - (C[0] - A[0]) * (B[2] - A[2])
        if abs(d) < 1e-12:
            continue
        l1 = ((X - A[0]) * (C[2] - A[2]) - (C[0] - A[0]) * (Z - A[2])) / d
        l2 = ((B[0] - A[0]) * (Z - A[2]) - (X - A[0]) * (B[2] - A[2])) / d
        ins = (l1 >= 0) & (l2 >= 0) & (l1 + l2 <= 1)
        yd = A[1] + l1 * (B[1] - A[1]) + l2 * (C[1] - A[1])
        upd = ins & (yd < best)                 # smaller y = nearer the viewer
        best[upd] = yd[upd]
        fid[upd] = i
    return fid, (x0, x0 + s * len(xs), z0, z0 + s * len(zs)), R


def draw_solid_exact(ax, V, F, normals, az, el, npx=900):
    fid, ext, R = raycast(V, F, normals, az, el, npx)
    rgb = np.zeros(fid.shape + (4,))
    for i in range(len(F)):
        c = matplotlib.colors.to_rgba(CAT[i])
        rgb[fid == i] = c
    edge = np.zeros(fid.shape, bool)               # outline wherever the visible face changes
    edge[:, 1:] |= fid[:, 1:] != fid[:, :-1]
    edge[1:, :] |= fid[1:, :] != fid[:-1, :]
    edge[:, :-1] |= fid[:, 1:] != fid[:, :-1]
    edge[:-1, :] |= fid[1:, :] != fid[:-1, :]
    from scipy import ndimage
    edge = ndimage.binary_dilation(edge, iterations=2) & (fid >= 0)
    rgb[edge] = matplotlib.colors.to_rgba(BG)
    ax.imshow(rgb, origin="lower", extent=ext, interpolation="nearest", zorder=2)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    ax.set_aspect("equal")
    ax.axis("off")
    return fid, ext, R


def enclosed_background(fid):
    """Background pixels not connected to the image border: places where one can see through the solid."""
    from scipy import ndimage
    lab, n = ndimage.label(fid < 0)
    border = set(np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
    inner = [k for k in range(1, n + 1) if k not in border]
    if not inner:
        return None
    k = max(inner, key=lambda k: (lab == k).sum())
    return np.argwhere(lab == k).mean(axis=0)          # (row, col)


def fig_hero(V, F, normals, reflex):
    convex = [i for i, r in enumerate(reflex) if not any(r)]
    assert len(convex) == 1
    fig = plt.figure(figsize=(12.5, 6.4))
    frame(fig)
    fig.text(0.03, 0.935, "Szilassi's polyhedron: seven faces, each touching all six others", fontsize=17,
             weight="bold")
    fig.text(0.03, 0.875, "A ring-shaped solid of seven flat hexagons (one colour each), drawn from the published exact "
             "coordinates. Three views.", fontsize=13, color=AXC)
    views = [(0, 20, "front, from slightly above"), (180, 20, "back, from slightly above"), (195, -35, "from below")]
    for j, (az, el, ttl) in enumerate(views):
        ax = fig.add_axes([0.025 + j * 0.245, 0.07, 0.235, 0.74])
        fid, ext, R = draw_solid_exact(ax, V, F, normals, az, el)
        ax.set_title(ttl, fontsize=13, color=AXC, loc="center")
        sx = (ext[1] - ext[0]) / fid.shape[1]
        sz = (ext[3] - ext[2]) / fid.shape[0]
        hole = enclosed_background(fid)
        if hole is not None and j < 2:
            p = (ext[0] + (hole[1] + 0.5) * sx, ext[2] + (hole[0] + 0.5) * sz)
            ax.annotate("see-through hole", p, xytext=(0, -150), textcoords="offset points", ha="center",
                        fontsize=12.5, color=TXT, arrowprops=dict(arrowstyle="-", color=TXT, lw=1.2, shrinkB=1),
                        bbox=dict(boxstyle="round,pad=0.25", fc=BG, ec="none"))
        if j == 2:
            rc = np.argwhere(fid == convex[0]).mean(axis=0)
            p = (ext[0] + (rc[1] + 0.5) * sx, ext[2] + (rc[0] + 0.5) * sz)
            ax.annotate("the one\nun-dented face", p, xytext=(72, -150), textcoords="offset points", ha="center",
                        fontsize=12.5, color=TXT, weight="bold", bbox=dict(boxstyle="round,pad=0.25", fc=BG, ec="none"),
                        arrowprops=dict(arrowstyle="-", color=TXT, lw=1.2, shrinkB=1))
    lx = fig.add_axes([0.78, 0.10, 0.20, 0.70])
    lx.axis("off")
    lx.set_xlim(0, 1)
    lx.set_ylim(-0.6, 7.6)
    lx.text(0, 7.3, "face: inward corners", fontsize=12.5, color=AXC, va="center")
    for i in range(7):
        y = 6.2 - i
        lx.add_patch(Polygon([(0.0, y - 0.3), (0.15, y - 0.3), (0.15, y + 0.3), (0.0, y + 0.3)], fc=CAT[i], ec="none"))
        k = sum(reflex[i])
        lx.text(0.22, y, "%d: %s" % (i + 1, "%d (dented)" % k if k else "0 (un-dented)"), fontsize=12.5,
                color=TXT, va="center", weight="bold" if not k else "normal")
    save(fig, "szilassi-polyhedron.svg")
    return convex


def fig_faces(V, F, normals, reflex):
    fig = plt.figure(figsize=(12.5, 5.0))
    frame(fig)
    fig.text(0.03, 0.91, "Six dented faces and exactly one un-dented face", fontsize=17, weight="bold")
    fig.text(0.03, 0.845, "The seven hexagons laid flat, all at the same scale. Open circles: corners that point inwards "
             "(%d of %d)." % (sum(map(sum, reflex)), 6 * len(F)), fontsize=13, color=AXC)
    shapes = []
    for i, f in enumerate(F):
        e1, e2, _ = plane_basis(normals[i])
        P = np.array([[np.array([float(x) for x in V[v]]) @ e1, np.array([float(x) for x in V[v]]) @ e2] for v in f])
        P -= P.mean(0)
        shapes.append(P)
    half = max(np.abs(np.vstack(shapes)).max(axis=0)) + 2
    for i, P in enumerate(shapes):
        ax = fig.add_axes([0.02 + i * 0.137, 0.06, 0.13, 0.66])
        ax.add_patch(Polygon(P, closed=True, fc=CAT[i], ec=TXT if not any(reflex[i]) else CAT[i], lw=2))
        for s in range(6):
            if reflex[i][s]:
                ax.plot(*P[s], "o", ms=10, mfc=BG, mec=TXT, mew=2, zorder=5)
        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)
        ax.set_aspect("equal")
        ax.axis("off")
        k = sum(reflex[i])
        ax.set_title("face %d\n%s" % (i + 1, "%d inward" % k if k else "un-dented"), fontsize=12.5,
                     color=TXT, weight="bold" if not k else "normal")
    save(fig, "szilassi-faces-flat.svg")


def fig_certs(root):
    rows = list(csv.DictReader(open(os.path.join(root, "results", "certificates.csv"), encoding="utf-8")))
    assert all(r["verdict"] == "VERIFIED" for r in rows)
    grp = lambda g: g.split("#")[0]
    names = {"C2": "half-turn", "C3": "one-third turn", "S4": "quarter-turn with mirror",
             "D2": "D2 (4 motions)", "D3": "D3 (6 motions)", "T": "T (12 motions)"}
    order = ["C2", "C3", "S4", "D2", "D3", "T"]
    fig = plt.figure(figsize=(12.5, 5.8))
    frame(fig)
    ax = fig.add_axes([0.25, 0.14, 0.63, 0.60])
    fig.text(0.03, 0.92, "Every possible symmetry ruled out: %d cases, %d checked proofs of impossibility"
             % (len(rows), sum(r["verdict"] == "VERIFIED" for r in rows)), fontsize=17, weight="bold")
    fig.text(0.03, 0.80, "One dot per case (a 12-face pattern with a given symmetry). Across: time the solver needed "
             "to prove it impossible.\nEvery proof file was re-checked by a separate checker program: all VERIFIED.", fontsize=13,
             color=AXC, linespacing=1.4)
    rng = np.random.default_rng(0)
    for y, g in enumerate(order):
        t = np.array([float(r["solve_s"]) for r in rows if grp(r["group"]) == g])
        jit = (rng.random(len(t)) - 0.5) * 0.36 if len(t) > 1 else np.zeros(1)
        ax.plot(t, y + jit, "o", ms=8, color=ACC, alpha=1.0, mec=BG, mew=0.8, zorder=3)
        ax.text(1.02, y, "%d case%s" % (len(t), "s" if len(t) > 1 else ""), transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=12.5, color=TXT, clip_on=False)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([names[g] for g in order])
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xscale("log")
    ticks = [5, 30, 60, 300, 1800, 3600, 7200]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["5 s", "30 s", "1 min", "5 min", "30 min", "1 h", "2 h"])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(4, 12000)
    ax.set_xlabel("solver time per case (log scale)")
    ax.tick_params(axis="y", length=0, labelsize=13, labelcolor=TXT)
    ax.grid(axis="y", visible=False)
    save(fig, "szilassi-symmetry-cases.svg")
    return rows


def main():
    root = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), ".."))
    style()
    V, F = read_table3(root)
    assert len(V) == 14 and len(F) == 7 and all(len(f) == 6 for f in F)
    normals, reflex, E = analyse(V, F)
    convex = fig_hero(V, F, normals, reflex)
    fig_faces(V, F, normals, reflex)
    rows = fig_certs(root)
    print("faces 7, shared edges %d, reflex corners %d of 42, convex faces %s (normal %s), certificates %d"
          % (len(E), sum(map(sum, reflex)), [c + 1 for c in convex], normals[convex[0]], len(rows)))
    print("inward corners per face:", [sum(r) for r in reflex])


if __name__ == "__main__":
    main()
