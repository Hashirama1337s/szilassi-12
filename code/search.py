"""Numerical realisation search (GPU, batched): minimise an exact-violation functional over plane coefficients.

For every test of checker A (S: face self-crossing, X: edge through a foreign face) we measure a signed depth:
  X(ij,k):  min(depth of v(i,j,k) inside segment ij, signed depth of v(i,j,k) inside face k)   > 0  <=> violated
  S(i,c,d): min(depth of v(i,c,d) inside segment ic, depth inside segment id)                   > 0  <=> violated
All lengths are divided by the configuration's scale (RMS radius of its vertices). Loss = sum relu(depth + margin)^2 plus
barriers against degenerate configurations (very short edges, nearly dependent normals at a vertex). A configuration
whose hard count (depth > 0) is zero is handed to the exact checker after rounding to rationals."""
import itertools, math, time, sys, os
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k12

try:
    import psutil
    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception:
    pass


class Problem:
    def __init__(self, S, dev='cuda'):
        self.S, self.dev = S, dev
        n = self.n = S.n
        self.triples = list(itertools.combinations(range(n), 3))
        tid = {}
        for q, t in enumerate(self.triples):
            for pm in itertools.permutations(t):
                tid[pm] = q
        self.tid = tid
        m = n - 1
        self.face_v = [[tid[(i, S.link[i][s], S.link[i][(s + 1) % m])] for s in range(m)] for i in range(n)]
        self.vert = sorted({tid[t] for t in S.tris})                      # the 44 (14) polyhedron vertices
        seg = {}
        for (i, j), (a, b) in S.elink.items():
            seg[(i, j)] = seg[(j, i)] = (tid[(i, j, a)], tid[(i, j, b)])
        self.seg = seg
        X = []
        for (i, j), (a, b) in S.elink.items():
            for k in range(n):
                if k not in (i, j, a, b):
                    X.append((tid[(i, j, k)], *seg[(i, j)], k))
        Sx = []
        for i in range(n):
            L = S.link[i]
            for s, t in itertools.combinations(range(m), 2):
                if (t - s) % m in (1, m - 1):
                    continue
                c, d = L[s], L[t]
                Sx.append((tid[(i, c, d)], *seg[(i, c)], *seg[(i, d)]))
        self.X = torch.tensor(X, device=dev)
        self.Sx = torch.tensor(Sx, device=dev)
        self.T = torch.tensor(self.triples, device=dev)
        self.FV = torch.tensor(self.face_v, device=dev)
        self.EV = torch.tensor(sorted({seg[(i, j)] for (i, j) in S.elink}), device=dev)
        self.VT = torch.tensor([self.triples[q] for q in self.vert], device=dev)
        self.vert_t = torch.tensor(self.vert, device=dev)

    # ---------------------------------------------------------------------------------------------------------
    def points(self, H):
        """H: (B, n, 4) -> unit-normal planes, all C(n,3) triple points (B, T, 3) and their normal determinants."""
        a = H[..., :3]
        nrm = a.norm(dim=-1, keepdim=True)
        a, d = a / nrm, H[..., 3] / nrm[..., 0]
        i, j, k = self.T[:, 0], self.T[:, 1], self.T[:, 2]
        ai, aj, ak = a[:, i], a[:, j], a[:, k]
        cjk, cki, cij = torch.cross(aj, ak, dim=-1), torch.cross(ak, ai, dim=-1), torch.cross(ai, aj, dim=-1)
        det = (ai * cjk).sum(-1)
        dets = torch.where(det.abs() < 1e-9, torch.full_like(det, 1e-9) * torch.sign(det + 1e-30), det)
        P = -(d[:, i, None] * cjk + d[:, j, None] * cki + d[:, k, None] * cij) / dets[..., None]
        return a, P, det

    @staticmethod
    def seg_depth(p, A, B):
        d = B - A
        L2 = (d * d).sum(-1).clamp_min(1e-18)
        s = ((p - A) * d).sum(-1) / L2
        return torch.minimum(s, 1 - s) * L2.sqrt()

    def face_depth(self, p, k, P, a):
        """signed depth of point p (B, N, 3), lying in plane k (N,), inside face k: +dist inside, -dist outside."""
        B = P.shape[0]
        fv = self.FV[k]                                                   # (N, m)
        Q = P[:, fv]                                                      # (B, N, m, 3)
        Q2 = torch.roll(Q, -1, dims=2)
        e = Q2 - Q
        w = p[:, :, None, :] - Q
        t = ((w * e).sum(-1) / (e * e).sum(-1).clamp_min(1e-18)).clamp(0, 1)
        dist = (w - t[..., None] * e).norm(dim=-1).min(-1).values        # (B, N)
        u, v = Q - p[:, :, None, :], Q2 - p[:, :, None, :]
        nk = a[:, k][:, :, None, :]
        ang = torch.atan2((torch.cross(u, v, dim=-1) * nk).sum(-1), (u * v).sum(-1)).sum(-1) / (2 * math.pi)
        inside = ang.abs() > 0.5
        return torch.where(inside, dist, -dist)

    def evaluate(self, H, margin=0.01, short=0.03, detmin=0.03):
        a, P, det = self.points(H)
        V = P[:, self.vert_t]
        scale = (V - V.mean(1, keepdim=True)).norm(dim=-1).pow(2).mean(1).sqrt().clamp_min(1e-9)   # (B,)
        X = self.X
        p, A, Bq = P[:, X[:, 0]], P[:, X[:, 1]], P[:, X[:, 2]]
        dX = torch.minimum(self.seg_depth(p, A, Bq), self.face_depth(p, X[:, 3], P, a)) / scale[:, None]
        Sx = self.Sx
        p = P[:, Sx[:, 0]]
        dS = torch.minimum(self.seg_depth(p, P[:, Sx[:, 1]], P[:, Sx[:, 2]]),
                           self.seg_depth(p, P[:, Sx[:, 3]], P[:, Sx[:, 4]])) / scale[:, None]
        el = (P[:, self.EV[:, 0]] - P[:, self.EV[:, 1]]).norm(dim=-1) / scale[:, None]
        vd = det[:, self.vert_t].abs()
        loss = (torch.relu(dX + margin) ** 2).sum(1) + (torch.relu(dS + margin) ** 2).sum(1) \
            + 10 * (torch.relu(short - el) ** 2).sum(1) + 10 * (torch.relu(detmin - vd) ** 2).sum(1)
        hard = (dX > 0).sum(1) + (dS > 0).sum(1)
        degen = (el < short / 3).any(1) | (vd < detmin / 3).any(1)
        return loss, hard, degen, dict(dX=dX, dS=dS)


def random_planes(B, n, dev, gen=None):
    a = torch.randn(B, n, 3, device=dev, generator=gen)
    a = a / a.norm(dim=-1, keepdim=True)
    d = -torch.ones(B, n, 1, device=dev) + 0.3 * torch.randn(B, n, 1, device=dev, generator=gen)
    return torch.cat([a, d], -1)


def run(S, B=2048, steps=3000, lr=0.02, restart_every=600, seed=0, time_limit=600, log=print, dev='cuda', noise=0.0,
        margin=0.01, short=0.03, detmin=0.03, param=None):
    """param: optional (init_fn, map_fn) for a reduced (e.g. symmetric) parametrisation: H = map_fn(Z)."""
    torch.manual_seed(seed)
    pb = Problem(S, dev)
    init_fn, map_fn = param if param else (lambda b: random_planes(b, S.n, dev), lambda Z: Z)
    H = init_fn(B).requires_grad_(True)
    opt = torch.optim.Adam([H], lr=lr)
    best = torch.full((B,), 1e9, device=dev)
    stall = torch.zeros(B, dtype=torch.long, device=dev)
    found, t0, it, starts = [], time.time(), 0, B
    minhard_seen = 10 ** 9
    near = []                                   # (hard, loss, H) snapshots of the closest misses
    while time.time() - t0 < time_limit:
        Hf = map_fn(H)
        loss, hard, degen, _ = pb.evaluate(Hf, margin, short, detmin)
        opt.zero_grad()
        loss.sum().backward()
        opt.step()
        with torch.no_grad():
            h = torch.where(degen, torch.full_like(hard, 10 ** 6), hard)
            hm = int(h.min())
            if hm <= max(minhard_seen, 2) and hm > 0 and len(near) < 200:
                q = int(torch.where(h == hm, loss, torch.full_like(loss, 1e30)).argmin())
                near.append((hm, float(loss[q]), Hf[q].detach().double().cpu().numpy().copy()))
            minhard_seen = min(minhard_seen, hm)
            ok = (h == 0).nonzero().flatten()
            for q in ok.tolist():
                found.append(Hf[q].detach().double().cpu().numpy().copy())
            if noise:
                H += noise * torch.randn_like(H) * (h > 0).view(-1, *([1] * (H.dim() - 1)))
            improved = loss < best * 0.999
            best = torch.where(improved, loss, best)
            stall = torch.where(improved, torch.zeros_like(stall), stall + 1)
            re = (stall > restart_every) | (h == 0)
            if re.any():
                idx = re.nonzero().flatten()
                H[idx] = init_fn(len(idx))
                best[idx] = 1e9
                stall[idx] = 0
                opt.state[H]['exp_avg'][idx] = 0
                opt.state[H]['exp_avg_sq'][idx] = 0
                starts += len(idx)
        it += 1
        if it % 500 == 0:
            log(f'  it {it:6d}  t {time.time() - t0:6.0f}s  starts {starts:8d}  min hard now {int(h.min()):3d}  '
                f'min hard seen {minhard_seen:3d}  found {len(found)}')
        if len(found) >= 5:
            break
    near.sort(key=lambda x: (x[0], x[1]))
    return found, dict(starts=starts, iters=it, minhard=minhard_seen, secs=time.time() - t0, near=near[:10])


if __name__ == '__main__':
    S7 = k12.Surface(k12.torus7(), 'torus7')
    found, info = run(S7, B=1024, time_limit=120)
    print(info, len(found))
