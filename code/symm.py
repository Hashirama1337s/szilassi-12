"""Reduced parametrisation: configurations invariant under a 180-degree rotation R about the z-axis that acts on the faces
by an involution perm (an orientation-preserving automorphism of the triangulation). Fixed faces get planes z = c."""
import torch

def c2_param(n, perm, dev='cuda'):
    fixed = [i for i in range(n) if perm[i] == i]
    reps = [i for i in range(n) if perm[i] > i]
    npar = 4 * len(reps) + len(fixed)
    R = torch.tensor([-1.0, -1.0, 1.0, 1.0], device=dev)
    def map_fn(Z):
        B = Z.shape[0]
        H = torch.zeros(B, n, 4, device=dev, dtype=Z.dtype)
        for q, i in enumerate(reps):
            h = Z[:, 4 * q:4 * q + 4]
            H[:, i] = h
            H[:, perm[i]] = h * R
        for q, i in enumerate(fixed):
            H[:, i, 2] = 1.0
            H[:, i, 3] = Z[:, 4 * len(reps) + q]
        return H
    def init_fn(B):
        Z = torch.randn(B, npar, device=dev)
        for q in range(len(reps)):
            a = Z[:, 4 * q:4 * q + 3]
            Z[:, 4 * q:4 * q + 3] = a / a.norm(dim=-1, keepdim=True)
            Z[:, 4 * q + 3] = -1 + 0.3 * torch.randn(B, device=dev)
        Z[:, 4 * len(reps):] = 0.5 * torch.randn(B, len(fixed), device=dev)
        return Z
    return init_fn, map_fn
