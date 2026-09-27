"""Symmetric dual SAT: python sym_sat.py KEY LABEL  (group from symgroup.candidate_groups; label T, D3, D2, S4, C3, C2[#k])"""
import sys, time, numpy as np
import k12, symgroup, run_sat
from pysat.solvers import Solver
try:
    import psutil; psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception: pass
key, lab = sys.argv[1], sys.argv[2]
extra = tuple(sys.argv[3].split(',')) if len(sys.argv) > 3 else ()
idx = int(lab.split('#')[1]) if '#' in lab else 0; lab = lab.split('#')[0]
S = {s.name.split('#')[-1]: s for s in k12.all59()}[key]
cands = [rho for l, rho in symgroup.candidate_groups(S.name) if l == lab]
rho = cands[idx]
group = [(g, int(round(np.linalg.det(M)))) for g, M in rho.items()]
t0 = time.time(); bd = run_sat.build('dual', S, group, extra); tb = time.time() - t0
with Solver(name='cadical195', bootstrap_with=bd.clauses) as s:
    r = s.solve()
print(f'sym #{key:3s} extra={",".join(extra) or "-"} {lab}{"#" + str(idx) if idx else "":3s} |G|={len(group):2d} {"SAT" if r else "UNSAT":6s} vars {bd.nv} clauses {len(bd.clauses)} build {tb:.1f}s solve {time.time() - t0 - tb:.1f}s', flush=True)
