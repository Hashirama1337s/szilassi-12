"""Exhaustive map of T-symmetric configurations of #54: plane a.x = 1 with unit normal a on a Fibonacci grid of the sphere."""
import math, torch, numpy as np, k12, symgroup, search
S = {s.name.split('_')[-1]: s for s in k12.all59()}['#54']
rho = [r for l, r in symgroup.candidate_groups(S.name) if l == 'T'][0]
init_fn, map_fn, reps = symgroup.param(rho)
pb = search.Problem(S)
N = 2_000_000
k = torch.arange(N, device='cuda', dtype=torch.float64) + 0.5
z = 1 - 2 * k / N; r = torch.sqrt(1 - z * z); ph = math.pi * (1 + 5 ** 0.5) * k
A = torch.stack([r * torch.cos(ph), r * torch.sin(ph), z], -1).float()
best = []; hist = torch.zeros(700, dtype=torch.long)
with torch.no_grad():
    for s in range(0, N, 8192):
        a = A[s:s + 8192]
        Z = torch.cat([a, -torch.ones(len(a), 1, device='cuda')], -1)[:, None, :]
        loss, hard, degen, _ = pb.evaluate(map_fn(Z), margin=0.0)
        h = torch.where(degen, torch.full_like(hard, 699), hard).clamp(max=699)
        hist += torch.bincount(h.cpu(), minlength=700)
        q = int(h.argmin()); best.append((int(h[q]), a[q].tolist()))
best.sort()
print('grid spacing ~', round(math.degrees(math.sqrt(4 * math.pi / N)), 4), 'deg; min hard', best[0][0], 'at', [round(x, 4) for x in best[0][1]])
print('histogram (hard: count) lowest:', [(i, int(c)) for i, c in enumerate(hist) if c > 0][:8])
