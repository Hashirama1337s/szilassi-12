import json, random, numpy as np, torch
from fractions import Fraction
import k12, search, check, export, symm
random.seed(7); torch.manual_seed(7)
S7 = k12.Surface(k12.torus7(), 'torus7')
L59 = k12.all59()
cases = []
f, _ = search.run(S7, B=1024, time_limit=20, param=symm.c2_param(7, [(-i) % 7 for i in range(7)]), log=lambda s: None)
cases += [(S7, h) for h in f[:3]]
f, info = search.run(S7, B=1024, time_limit=40, noise=0.003, log=lambda s: None, seed=3)
cases += [(S7, h) for h in f[:2]] + [(S7, x[2]) for x in info['near'][:4]]
for k in (0, 53):
    f, info = search.run(L59[k], B=1024, time_limit=30, noise=0.003, log=lambda s: None, seed=5)
    cases += [(L59[k], x[2]) for x in info['near'][:3]]
cases += [(S7, search.random_planes(1, 7, 'cpu')[0].numpy())]
random.shuffle(cases)
key = {}
for q, (S, H) in enumerate(cases):
    Hq = [[Fraction(round(float(x) * 10**6), 10**6) for x in h] for h in H]
    try: ok = check.check(S, Hq)[0]
    except ValueError as e: ok = 'degenerate'
    json.dump(export.polyhedron(S, Hq), open(f'checkerB_tests/case{q:02d}.json', 'w'))
    key[f'case{q:02d}'] = ok
json.dump(key, open('checkerB_tests/expected_checkerA.json', 'w'))
print({k: v for k, v in key.items()})
