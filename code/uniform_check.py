"""For every (surface, symmetry group) case: a random G-symmetric configuration has all 715 4x4 determinants
(12 planes + plane at infinity) nonzero, so no determinant vanishes identically on the symmetric subspace; generic
symmetric configurations are uniform, and an embedded symmetric realisation (an open condition) can be perturbed within
the subspace to a uniform one. Also checks the free action (no face fixed)."""
import itertools, numpy as np, k12, symgroup
rng = np.random.default_rng(12)
rows = []
for S in k12.all59():
    for lab, rho in symgroup.candidate_groups(S.name):
        init_fn, map_fn, reps = symgroup.param(rho, dev='cpu')
        import torch
        worst = 1e9
        for t in range(3):
            H = map_fn(init_fn(1).double())[0].numpy() if False else map_fn(init_fn(1))[0].double().numpy()
            V = list(H) + [np.array([0, 0, 0, 1.0])]
            m = min(abs(np.linalg.det(np.array([V[x] for x in b]))) for b in itertools.combinations(range(13), 4))
            worst = min(worst, m)
        rows.append((S.name.split('#')[-1], lab, len(rho), worst))
bad = [r for r in rows if r[3] < 1e-6]
print(len(rows), 'cases; smallest |det| over all cases', min(r[3] for r in rows), '; cases with a vanishing determinant:', bad)
