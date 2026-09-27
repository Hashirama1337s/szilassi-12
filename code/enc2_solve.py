"""enc2_solve -- build and solve enc2 instances.
  python enc2_solve.py one <cnf> <out.json>            solve one DIMACS file with CaDiCaL 1.9.5 (pysat) -> json
  python enc2_solve.py t7                              the 7-face torus instance (must be SAT)
  python enc2_solve.py sym [key:label ...]             symmetric 12-face instances, 2 parallel, 40-min cap each
Solves run in child processes at below-normal priority; the cap is enforced by killing the child."""
import sys, os, time, json, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CAP = 2400
TARGETS = [('#54', 'T'), ('#3', 'D3'), ('#48', 'D3'), ('#59', 'D2'), ('#5', 'S4'), ('#45', 'C3')]


def write_cnf(E, path):
    cls = E.cnf()
    with open(path, 'w') as fh:
        fh.write('p cnf %d %d\n' % (E.nvar, len(cls)))
        for c in cls:
            fh.write(' '.join(map(str, c)) + ' 0\n')
    return len(cls)


def solve_one(cnf, out):
    from pysat.formula import CNF
    from pysat.solvers import Solver
    F = CNF(from_file=cnf)
    t0 = time.time()
    with Solver(name='cadical195', bootstrap_with=F.clauses) as s:
        r = s.solve()
        model = s.get_model() if r else None
        st = s.accum_stats()
    json.dump(dict(status='SAT' if r else 'UNSAT', time=time.time() - t0, model=model, stats=st), open(out, 'w'))


def launch(cnf, out):
    flags = subprocess.BELOW_NORMAL_PRIORITY_CLASS if os.name == 'nt' else 0
    return subprocess.Popen([sys.executable, os.path.abspath(__file__), 'one', cnf, out], creationflags=flags,
                            stdout=subprocess.DEVNULL, stderr=open(out + '.err', 'w'))


def check_model(E, model):
    """the model's chi values satisfy every condition (evaluated through the gate definitions)."""
    val = np.zeros(E.nchi, dtype=bool)
    for l in model:
        if 1 <= abs(l) <= E.nchi:
            val[abs(l) - 1] = l > 0
    V = E.evaluate(val[:, None])
    idx, viol = E.violated(V, E.fams)
    return int(viol.sum())


def sign_fixes(E, orbits):
    """Symmetry breaking (justification): negating h_i on a whole G-orbit, or the central inversion x -> -x (which
    commutes with every M and negates chi), maps G-symmetric embedded configurations to G-symmetric embedded ones and
    changes chi(B) by prod_orbits (-1)^{|B n O|} resp. by -1; choosing k+1 bases whose GF(2) vectors
    (|B n O_1|,...,|B n O_k|, 1) are independent, we may therefore assume these k+1 bases are positive."""
    import itertools
    k = len(orbits)
    chosen, rows = [], []

    def rank(M):
        M = [r[:] for r in M]
        rk = 0
        for c in range(k + 1):
            piv = next((r for r in range(rk, len(M)) if M[r][c]), None)
            if piv is None:
                continue
            M[rk], M[piv] = M[piv], M[rk]
            for r in range(len(M)):
                if r != rk and M[r][c]:
                    M[r] = [x ^ y for x, y in zip(M[r], M[rk])]
            rk += 1
        return rk
    for B in E.bases:
        vec = [len(set(B) & O) % 2 for O in orbits] + [1]
        if rank(rows + [vec]) > len(rows):
            rows.append(vec)
            chosen.append(B)
        if len(rows) == k + 1:
            break
    assert len(rows) == k + 1
    return [[E.chi(B)] for B in chosen], chosen


def build(key, label, sb=False):
    import symgroup
    from enc2_core import Enc, surface
    S = surface(key)
    rho = dict(symgroup.candidate_groups(S.name))[label]
    sym = [(g, int(round(np.linalg.det(M)))) for g, M in rho.items()]
    E = Enc(S, sym=sym)
    path = os.path.join(HERE, 'enc2_cnf_%s_%s%s.cnf' % (key.strip('#'), label, '_sb' if sb else ''))
    if sb:
        orbits, seen = [], set()
        for f in range(S.n):
            if f not in seen:
                O = {g[f] for g in rho}
                seen |= O
                orbits.append(O)
        units, chosen = sign_fixes(E, orbits)
        E.defs += units
        print('   symmetry breaking: %d orbits, fixed positive bases %s' % (len(orbits), chosen), flush=True)
    nc = write_cnf(E, path)
    nrep = len({abs(E.rep[v]) for v in range(1, E.nchi + 1)})
    return E, path, nc, nrep


def run_sym(targets, sb=False):
    jobs = []
    for key, label in targets:
        E, path, nc, nrep = build(key, label, sb=sb)
        print('%s %s: chi orbits %d, vars %d, clauses %d, trivial contradiction %s' % (
            key, label, nrep, E.nvar, nc, E.contradiction), flush=True)
        jobs.append((key, label, E, path))
    pending, running, results = list(jobs), [], {}
    while pending or running:
        while pending and len(running) < 2:
            key, label, E, path = pending.pop(0)
            out = path[:-4] + '.json'
            if os.path.exists(out):
                os.remove(out)
            running.append((key, label, E, out, launch(path, out), time.time()))
        time.sleep(2)
        for job in list(running):
            key, label, E, out, p, t0 = job
            el = time.time() - t0
            if p.poll() is not None:
                running.remove(job)
                R = json.load(open(out)) if os.path.exists(out) else dict(status='ERROR', time=el)
                extra = ''
                if R['status'] == 'SAT':
                    extra = ' model check: %d violated conditions' % check_model(E, R['model'])
                results[(key, label)] = (R['status'], R['time'])
                print('RESULT %s %s: %s in %.1f s%s' % (key, label, R['status'], R['time'], extra), flush=True)
            elif el > CAP:
                p.kill()
                running.remove(job)
                results[(key, label)] = ('TIMEOUT', el)
                print('RESULT %s %s: TIMEOUT after %.0f s' % (key, label, el), flush=True)
    return results


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'one':
        solve_one(sys.argv[2], sys.argv[3])
    elif cmd == 't7':
        from enc2_core import Enc, surface
        E = Enc(surface('t7'))
        path = os.path.join(HERE, 'enc2_cnf_t7.cnf')
        nc = write_cnf(E, path)
        out = os.path.join(HERE, 'enc2_cnf_t7.json')
        p = launch(path, out)
        p.wait()
        R = json.load(open(out))
        print('t7: vars %d clauses %d -> %s in %.2f s; model check: %d violated' % (
            E.nvar, nc, R['status'], R['time'], check_model(E, R['model']) if R['model'] else -1))
    elif cmd in ('sym', 'symsb'):
        tg = [tuple(a.split(':')) for a in sys.argv[2:]] or TARGETS
        run_sym(tg, sb=(cmd == 'symsb'))
