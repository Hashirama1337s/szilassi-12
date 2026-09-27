"""Custom-weighted GPU search for near-miss 12-face Szilassi configurations: heavily penalise genuine edge-face
piercings (X), lightly penalise face self-crossings (S) -- the goal is the best few-crossing near-miss, not a
full embedding. Reuses search.Problem for the (validated) geometry; never edits search.py or check.py."""
import itertools, time, os, sys
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import search as search_mod
import k12


def face_idx_of_Sx(S):
    """Same order as search.Problem.__init__'s Sx list; index i of the face each S-test belongs to."""
    n = S.n
    m = n - 1
    out = []
    for i in range(n):
        for s, t in itertools.combinations(range(m), 2):
            if (t - s) % m in (1, m - 1):
                continue
            out.append(i)
    return out


def make_problem(S, dev='cuda'):
    pb = search_mod.Problem(S, dev)
    fidx = face_idx_of_Sx(S)
    assert len(fidx) == pb.Sx.shape[0]
    pb.face_idx = torch.tensor(fidx, device=dev, dtype=torch.long)
    return pb


def evaluate_weighted(pb, H, wX=40.0, wS=1.0, wConc=0.0, margin=0.01, short=0.03, detmin=0.03, face_weight=None):
    """face_weight: optional (n,) tensor of per-face multipliers on the S loss (gathered via pb.face_idx), used
    to bias the optimiser to concentrate self-crossings onto a chosen subset of faces."""
    a, P, det = pb.points(H)
    V = P[:, pb.vert_t]
    scale = (V - V.mean(1, keepdim=True)).norm(dim=-1).pow(2).mean(1).sqrt().clamp_min(1e-9)
    X = pb.X
    p, A, Bq = P[:, X[:, 0]], P[:, X[:, 1]], P[:, X[:, 2]]
    dX = torch.minimum(pb.seg_depth(p, A, Bq), pb.face_depth(p, X[:, 3], P, a)) / scale[:, None]
    Sx = pb.Sx
    p2 = P[:, Sx[:, 0]]
    dS = torch.minimum(pb.seg_depth(p2, P[:, Sx[:, 1]], P[:, Sx[:, 2]]),
                        pb.seg_depth(p2, P[:, Sx[:, 3]], P[:, Sx[:, 4]])) / scale[:, None]
    el = (P[:, pb.EV[:, 0]] - P[:, pb.EV[:, 1]]).norm(dim=-1) / scale[:, None]
    vd = det[:, pb.vert_t].abs()
    sw = wS if face_weight is None else wS * face_weight[pb.face_idx][None, :]
    loss = wX * (torch.relu(dX + margin) ** 2).sum(1) + (sw * torch.relu(dS + margin) ** 2).sum(1) \
        + 10 * (torch.relu(short - el) ** 2).sum(1) + 10 * (torch.relu(detmin - vd) ** 2).sum(1)
    hardX = (dX > 0).sum(1)
    hardS = (dS > 0).sum(1)
    violated = (dS > 0).float()
    Bn = H.shape[0]
    face_max = torch.zeros(Bn, pb.S.n, device=H.device)
    face_max.scatter_reduce_(1, pb.face_idx.unsqueeze(0).expand(Bn, -1), violated, reduce='amax')
    facesS = (face_max > 0).sum(1)
    if wConc:
        # soft surrogate: encourage violations to concentrate into few faces (sum over faces of soft-OR of
        # sigmoid(depth)); only used in the refine stage, light weight.
        soft = torch.sigmoid(dS * 50.0)
        face_soft = torch.zeros(Bn, pb.S.n, device=H.device)
        face_soft.scatter_reduce_(1, pb.face_idx.unsqueeze(0).expand(Bn, -1), soft, reduce='amax')
        loss = loss + wConc * face_soft.sum(1)
    degen = (el < short / 3).any(1) | (vd < detmin / 3).any(1)
    return loss, hardX, hardS, facesS, degen


def random_init(S, dev):
    return lambda b: search_mod.random_planes(b, S.n, dev)


def seeded_init(base_H, dev='cuda', noise_scale=0.05, jitter_scale=0.0, noise_range=None):
    """base_H: (12,4) or (K,12,4) numpy array (a near-miss snapshot, possibly symmetric). Breaks symmetry by
    adding independent per-plane gaussian noise; jitter_scale additionally rescales each plane a little.
    noise_range=(lo,hi): if given, each of the B samples draws its own noise scale uniformly in [lo,hi]
    (log-spaced) instead of the fixed noise_scale, for diversity in one batch."""
    base = torch.tensor(np.asarray(base_H), dtype=torch.float32, device=dev)
    if base.dim() == 2:
        base = base.unsqueeze(0)
    K = base.shape[0]

    def init_fn(B):
        idx = torch.randint(0, K, (B,), device=dev)
        H = base[idx].clone()
        if noise_range:
            lo, hi = noise_range
            u = torch.rand(B, 1, 1, device=dev)
            scale = torch.exp(u * (np.log(hi) - np.log(lo)) + np.log(lo))
        else:
            scale = noise_scale
        H = H + scale * torch.randn_like(H)
        if jitter_scale:
            H = H * (1 + jitter_scale * torch.randn(H.shape[0], H.shape[1], 1, device=dev))
        return H
    return init_fn


def run(S, B=4096, lr=0.02, restart_every=500, seed=0, time_limit=60, dev='cuda',
        wX=40.0, wS=1.0, wConc=0.0, margin=0.01, short=0.03, detmin=0.03, init_fn=None, log=None,
        keep_top=8, face_weight=None, noise=0.0, noise_on='hardX', map_fn=None):
    """noise: if >0, continual gaussian kicks added each step to samples still failing `noise_on`
    ('hardX': only piercing-violated samples; 'any': hardX>0 or hardS>0), for basin-hopping past local minima.
    map_fn: optional Z -> H map (e.g. symgroup.param's map_fn) for searching a REDUCED (symmetric) parametrisation
    Z of shape (B, m, 4); init_fn must then produce Z, not the full (B, n, 4) H. Snapshots always store the full,
    mapped H (12 planes) since that is what check.py / kp_common understand."""
    torch.manual_seed(seed)
    pb = make_problem(S, dev)
    if init_fn is None:
        init_fn = random_init(S, dev)
    if map_fn is None:
        map_fn = lambda z: z
    Z = init_fn(B).requires_grad_(True)
    opt = torch.optim.Adam([Z], lr=lr)
    best = torch.full((B,), 1e18, device=dev)
    stall = torch.zeros(B, dtype=torch.long, device=dev)
    t0 = time.time()
    it = 0
    starts = B
    fw = torch.tensor(face_weight, dtype=torch.float32, device=dev) if face_weight is not None else None
    top = []  # list of (key, snapshot) kept sorted, small pool of diverse near-misses
    while time.time() - t0 < time_limit:
        H = map_fn(Z)
        loss, hardX, hardS, facesS, degen = evaluate_weighted(pb, H, wX, wS, wConc, margin, short, detmin, fw)
        opt.zero_grad()
        loss.sum().backward()
        opt.step()
        with torch.no_grad():
            keyval = hardX.double() * 1e6 + facesS.double() * 1e3 + hardS.double()
            keyval = torch.where(degen, torch.full_like(keyval, 1e15), keyval)
            qmin = int(torch.argmin(keyval))
            kmin = (int(hardX[qmin]), int(facesS[qmin]), int(hardS[qmin]))
            if (not top) or kmin < top[0][0]:
                snap = H[qmin].detach().double().cpu().numpy().copy()
                top.append((kmin, snap))
                top.sort(key=lambda x: x[0])
                top = top[:keep_top]
            if noise:
                bad = (hardX > 0) if noise_on == 'hardX' else ((hardX > 0) | (hardS > 0))
                Z += noise * torch.randn_like(Z) * bad.view(-1, *([1] * (Z.dim() - 1)))
            improved = loss < best * 0.999
            best = torch.where(improved, loss, best)
            stall = torch.where(improved, torch.zeros_like(stall), stall + 1)
            re = stall > restart_every
            if re.any():
                idx = re.nonzero().flatten()
                Z[idx] = init_fn(len(idx))
                best[idx] = 1e18
                stall[idx] = 0
                opt.state[Z]['exp_avg'][idx] = 0
                opt.state[Z]['exp_avg_sq'][idx] = 0
                starts += len(idx)
        it += 1
        if log and it % 1000 == 0:
            log(f'  it {it} t {time.time()-t0:.0f}s best_key {top[0][0] if top else None} starts {starts}')
    with torch.no_grad():
        # also grab a DIVERSE slice of the FINAL batch (not just the single historical-best trace), so callers
        # doing staged/annealed search have several distinct basins to seed the next phase from.
        H = map_fn(Z)
        loss, hardX, hardS, facesS, degen = evaluate_weighted(pb, H, wX, wS, wConc, margin, short, detmin, fw)
        keyval = hardX.double() * 1e6 + facesS.double() * 1e3 + hardS.double()
        keyval = torch.where(degen, torch.full_like(keyval, 1e15), keyval)
        order = torch.argsort(keyval)[:max(keep_top, 30)]
        pool = [((int(hardX[q]), int(facesS[q]), int(hardS[q])), H[q].detach().double().cpu().numpy().copy())
                for q in order.tolist()]
    return top, dict(iters=it, starts=starts, secs=time.time() - t0, pool=pool)


def symmetric_search_setup(S, label, dev='cuda'):
    """Return (init_fn, map_fn, reps, rho) for symmetry group `label` (e.g. 'C2','C3','S4',...) on surface S,
    via symgroup.candidate_groups; None if that group is not available for S."""
    import symgroup
    for lab, rho in symgroup.candidate_groups(S.name, n=S.n):
        if lab == label:
            init_fn, map_fn, reps = symgroup.param(rho, n=S.n, dev=dev)
            return init_fn, map_fn, reps, rho
    return None


if __name__ == '__main__':
    import k12
    S = k12.all59()[4]  # #5
    top, info = run(S, B=2048, time_limit=30, log=print)
    print(info)
    for k, snap in top[:5]:
        print(k)
