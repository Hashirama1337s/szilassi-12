import search, k12, symm, check, numpy as np, sys
from fractions import Fraction
S7 = k12.Surface(k12.torus7(), 'torus7')
mode = sys.argv[1]
if mode == 'noise':
    found, info = search.run(S7, B=2048, time_limit=90, noise=0.003, log=lambda s: None)
else:
    perm = [(-i) % 7 for i in range(7)]
    found, info = search.run(S7, B=2048, time_limit=90, param=symm.c2_param(7, perm), log=lambda s: None)
print(mode, info, 'found', len(found))
for H in found[:3]:
    Hq = [[Fraction(round(x * 10**6), 10**6) for x in h] for h in H]
    print(check.check(S7, Hq)[0], np.round(H, 4).tolist())
