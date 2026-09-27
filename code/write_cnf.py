"""Write a DIMACS CNF for a dual instance. usage: python write_cnf.py KEY LABEL [EXTRA]  (LABEL 'full' or a group label)"""
import sys, os, numpy as np
import k12, symgroup, run_sat
key, lab = sys.argv[1], sys.argv[2]
extra = tuple(sys.argv[3].split(',')) if len(sys.argv) > 3 and sys.argv[3] != '-' else ()
S = {s.name.split('#')[-1]: s for s in k12.all59()}[key]
group = None
if lab != 'full':
    idx = int(lab.split('#')[1]) if '#' in lab else 0
    rho = [r for l, r in symgroup.candidate_groups(S.name) if l == lab.split('#')[0]][idx]
    group = [(g, int(round(np.linalg.det(M)))) for g, M in rho.items()]
bd = run_sat.build('dual', S, group, extra)
os.makedirs('cnf', exist_ok=True)
fn = f'cnf/dual_{key}_{lab.replace("#", "-")}{"_" + "+".join(extra) if extra else ""}.cnf'
with open(fn, 'w') as f:
    f.write(f'p cnf {bd.nv} {len(bd.clauses)}\n')
    f.write(''.join(' '.join(map(str, c)) + ' 0\n' for c in bd.clauses))
print(fn)
