"""Which S/X constraints does a certified UNSAT proof actually use? usage: python core_analysis.py KEY LABEL [EXTRA]"""
import sys, os, subprocess, numpy as np, collections
import k12, symgroup, run_sat
key, lab = sys.argv[1], sys.argv[2]
extra = tuple(sys.argv[3].split(',')) if len(sys.argv) > 3 else ()
S = {s.name.split('#')[-1]: s for s in k12.all59()}[key]
idx = int(lab.split('#')[1]) if '#' in lab else 0
rho = [r for l, r in symgroup.candidate_groups(S.name) if l == lab.split('#')[0]][idx]
group = [(g, int(round(np.linalg.det(M)))) for g, M in rho.items()]
bd = run_sat.build('dual', S, group, extra)
cnf = f'cnf/core_{key}_{lab}.cnf'; prf = f'proofs/core_{key}_{lab}.drat'; core = f'proofs/core_{key}_{lab}.core'
with open(cnf, 'w') as f:
    f.write(f'p cnf {bd.nv} {len(bd.clauses)}\n' + ''.join(' '.join(map(str, c)) + ' 0\n' for c in bd.clauses))
subprocess.run([r'tools\cadical195\cadical.exe', '-q', cnf, prf], capture_output=True)
v = subprocess.run([r'tools\drat-trim\drat-trim.exe', cnf, prf, '-c', core, '-t', '40000'], capture_output=True, text=True)
print([l for l in v.stdout.splitlines() if l.startswith('s ') or 'core' in l][:3])
coreset = set()
for line in open(core):
    if line.startswith('p') or not line.strip(): continue
    coreset.add(tuple(sorted(int(x) for x in line.split()[:-1])))
marks = bd.marks
forbid = {}
for (st, label), (nx, _) in zip(marks, marks[1:]):
    forbid[tuple(sorted(bd.clauses[nx - 1]))] = label           # the constraint's own clause is the last in its range
used = [lab_ for cl, lab_ in forbid.items() if cl in coreset]
cnt = collections.Counter(u[0] for u in used)
print(f'#{key} {lab}: constraints used {dict(cnt)} of S {sum(1 for l in forbid.values() if l[0]=="S")}, X {sum(1 for l in forbid.values() if l[0]=="X")}')
faces = collections.Counter(f for u in used for f in u[1:])
print('faces involved:', sorted(faces.items()))
print('used:', sorted(used)[:40])
for p in (prf, core): os.remove(p)
