"""Build and solve the OM SAT instances.  usage: python run_sat.py MODE SURFACES SECONDS WORKERS TAG
MODE = dual | vertex ; SURFACES = t7 | all | comma list of Lutz numbers (e.g. 1,54).  Results appended to sat_<TAG>.log;
satisfying chirotopes saved to chi_<mode>_<surface>.json; every CNF is written to cnf/<mode>_<surface>.cnf for proof runs."""
import sys, os, time, json, threading, itertools
import multiprocessing as mp
import k12, omsat

try:
    import psutil
    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception:
    pass


def build(mode, S, group=None, extra=()):
    """group (dual only): list of (perm, s) with s = det of the isometry; imposes chi(g B) = s chi(B), INF fixed.
    Reorientation normalisation is then omitted (it could clash with the symmetric representative)."""
    n = S.n
    if mode == 'dual':
        bd = omsat.Builder(n + 1, 'cnf')
        omsat.gp3(bd, list(range(n + 1)))
        omsat.dual_constraints(bd, S, n, extra)
        if group is None:
            omsat.reorientation_fix(bd, list(range(n + 1)))
        else:
            for perm, s in group:
                g = list(perm) + [n]
                for B in bd.bases:
                    x, y = bd.neg(*B), bd.neg(*[g[e] for e in B])       # y must equal x xor (s == -1)
                    if s < 0:
                        y = -y
                    if set(g[e] for e in B) == set(B) and x != y and abs(x) == abs(y):
                        raise ValueError('symmetry forces a zero determinant (non-uniform)')
                    bd.clauses += [[-x, y], [x, -y]]
            if 'orbfix' in extra:
                orbit_reorientation_fix(bd, group, n)
    else:
        bd = omsat.Builder(n, 'cnf')
        omsat.gp3(bd, list(range(n)))
        omsat.acyclic(bd, list(range(n)))
        omsat.vertex_constraints(bd, S)
    return bd


def orbit_reorientation_fix(bd, group, n):
    """Symmetric instances: reorienting a whole face orbit (h -> -h on every plane of the orbit) or INF keeps the
    configuration G-invariant and every constraint unchanged, so these flips are symmetries of the instance. Fix
    chi = + on bases whose flip-parity vectors are GF(2)-independent (one representative per flip class)."""
    orbits = []
    seen = set()
    for f in range(n):
        if f not in seen:
            orb = sorted({g[f] for g, s in group})
            seen |= set(orb)
            orbits.append(orb)
    orbits.append([n])                                   # INF
    unit = {e: q for q, orb in enumerate(orbits) for e in orb}
    piv, picked = {}, []
    for B in bd.bases:
        w = 0
        for e in B:
            w ^= 1 << unit[e]
        while w:
            p = w.bit_length() - 1
            if p in piv:
                w ^= piv[p]
            else:
                piv[p] = w; picked.append(B); break
    for B in picked:
        bd.clauses.append([bd.var[B]])
    return picked


def job(args):
    mode, key, secs = args
    from pysat.solvers import Solver
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass
    S = k12.Surface(k12.torus7(), 't7') if key == 't7' else {s.name.split('#')[-1]: s for s in k12.all59()}[key]
    t0 = time.time()
    bd = build(mode, S)
    os.makedirs('cnf', exist_ok=True)
    with open(f'cnf/{mode}_{key}.cnf', 'w') as f:
        f.write(f'p cnf {bd.nv} {len(bd.clauses)}\n')
        for c in bd.clauses:
            f.write(' '.join(map(str, c)) + ' 0\n')
    tb = time.time() - t0
    with Solver(name='cadical195', bootstrap_with=bd.clauses) as s:
        timer = threading.Timer(secs, s.interrupt)
        timer.start()
        r = s.solve_limited(expect_interrupt=True)
        timer.cancel()
        model = s.get_model() if r else None
        st = s.accum_stats()
    res = {True: 'SAT', False: 'UNSAT', None: 'UNKNOWN'}[r]
    if r:
        m = set(x for x in model if x > 0)
        chi = {','.join(map(str, b)): (1 if bd.var[b] in m else -1) for b in bd.bases}
        json.dump(chi, open(f'chi_{mode}_{key}.json', 'w'))
    return mode, key, res, bd.nv, len(bd.clauses), tb, time.time() - t0 - tb, st.get('conflicts', 0)


if __name__ == '__main__':
    mode, surfs, secs, workers, tag = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    keys = ['t7'] if surfs == 't7' else ([str(q) for q in range(1, 60)] if surfs == 'all' else surfs.split(','))
    log = open(f'sat_{tag}.log', 'a', encoding='utf-8')
    log.write(f'=== {time.ctime()} {mode} {surfs} {secs}s x{workers}\n'); log.flush()
    with mp.Pool(workers) as pool:
        for mode_, key, res, nv, nc, tb, ts, conf in pool.imap_unordered(job, [(mode, k, secs) for k in keys]):
            line = f'{mode_:6s} #{key:3s} {res:7s} vars {nv:7d} clauses {nc:8d} build {tb:5.1f}s solve {ts:8.1f}s conflicts {conf}'
            print(line, flush=True); log.write(line + '\n'); log.flush()
    log.write('DONE\n'); log.close()
