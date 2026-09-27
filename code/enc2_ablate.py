"""enc2_ablate -- sanity controls for the symmetric instances: GP + symmetry alone must be SAT (symmetric uniform
chirotopes exist), and which condition families are needed for UNSAT on the fast #54/T instance."""
import sys, os, time
import numpy as np
import psutil

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from enc2_core import Enc, surface
from pysat.solvers import Solver
import symgroup


def solve(E, cap_note=''):
    t0 = time.time()
    with Solver(name='cadical195', bootstrap_with=E.cnf()) as s:
        r = s.solve()
    return ('SAT' if r else 'UNSAT'), time.time() - t0


if __name__ == '__main__':
    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    targets = [('#54', 'T'), ('#3', 'D3'), ('#48', 'D3'), ('#59', 'D2'), ('#5', 'S4'), ('#45', 'C3')]
    for key, lab in targets:
        S = surface(key)
        rho = dict(symgroup.candidate_groups(S.name))[lab]
        sym = [(g, int(round(np.linalg.det(M)))) for g, M in rho.items()]
        E = Enc(S, fams=('GP',), sym=sym)
        print(key, lab, 'GP+symmetry only:', *solve(E), flush=True)
    S = surface('#54')
    rho = dict(symgroup.candidate_groups(S.name))['T']
    sym = [(g, int(round(np.linalg.det(M)))) for g, M in rho.items()]
    for fams in [('GP', 'SIMP'), ('GP', 'PIERCE'), ('GP', 'REFL', 'LOC'), ('GP', 'SIMP', 'PIERCE'),
                 ('GP', 'SIMP', 'REFL', 'PIERCE'), ('SIMP', 'REFL', 'PIERCE', 'LOC')]:
        E = Enc(S, fams=fams, sym=sym)
        print('#54 T', '+'.join(fams), *solve(E), flush=True)
