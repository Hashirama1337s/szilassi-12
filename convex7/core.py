"""Core of the 7-face convexity proof (functions only; run convex7.py).
Every 7-face hexagonal polyhedral torus without overarching faces has exactly one convex face."""
import sys, os, itertools, json, hashlib, subprocess, time, random
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'code'))
import numpy as np
import k12, omsat, check, enc2_core

S7 = k12.Surface(k12.torus7(), 'torus7')
N, INF = 7, 7
CAD = os.environ.get('CADICAL', 'cadical')
DT = os.environ.get('DRAT_TRIM', 'drat-trim')
WORK = os.path.join(HERE, 'out')
os.makedirs(WORK, exist_ok=True)
OUT = os.path.join(WORK, 'results.jsonl')

# Grunbaum & Szilassi 2009, Appendix Table 3 (the Szilassi torus of Figure 10)
T3V = {1: (-24, 0, 24), 2: (24, 0, 24), 3: (0, Q('-25.2'), -24), 4: (0, Q('25.2'), -24), 5: (4, -10, -16),
       6: (-4, 10, -16), 7: (Q('-7.5'), Q('-7.5'), -6), 8: (Q('7.5'), Q('7.5'), -6), 9: (9, -5, 4), 10: (-9, 5, 4),
       11: (-14, 0, 4), 12: (14, 0, 4), 13: (-14, -5, 4), 14: (14, 5, 4)}
T3F = [(1, 2, 14, 10, 8, 6), (1, 6, 4, 3, 11, 13), (3, 5, 7, 8, 10, 11), (4, 6, 8, 7, 9, 12), (2, 5, 3, 4, 12, 14),
       (13, 9, 7, 5, 2, 1), (9, 13, 11, 10, 14, 12)]


def log(rec):
    rec['time'] = time.strftime('%H:%M:%S')
    with open(OUT, 'a', encoding='utf-8') as f:
        f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec))


def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def cross(a, b): return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
def dot(a, b): return sum(x * y for x, y in zip(a, b))


# ---------------------------------------------------------------- C0: the 7-face map is unique
def c0():
    V = range(7)
    tris = list(itertools.combinations(V, 3))
    edges = list(itertools.combinations(V, 2))
    found = []

    def rec(chosen, cnt):
        # every edge must lie in exactly 2 triangles (closed surface with all 21 edges: V-E+F = 7-21+14 = 0)
        e = next((e for e in edges if cnt[e] < 2), None)
        if e is None:
            found.append(list(chosen)); return
        cand = [t for t in tris if t not in chosen and e[0] in t and e[1] in t]
        # add the missing triangles on edge e all at once (a pair when e is still bare), so no set is found twice
        for add in itertools.combinations(cand, 2 - cnt[e]):
            es = [x for t in add for x in itertools.combinations(t, 2)]
            if any(cnt[x] + es.count(x) > 2 for x in es):
                continue
            for x in es: cnt[x] += 1
            chosen.extend(add)
            rec(chosen, cnt)
            del chosen[len(chosen) - len(add):]
            for x in es: cnt[x] -= 1

    rec([], {e: 0 for e in edges})
    sets = {tuple(sorted(s)) for s in found}
    good = []
    for s in sets:                       # vertex links must be single cycles (a 2-manifold)
        try:
            Sx = k12.Surface(list(s), 'x')
            good.append((s, Sx.genus))
        except (AssertionError, KeyError):
            pass
    ref = set(S7.tris)
    iso = 0
    for s, g in good:
        ss = set(s)
        if any({tuple(sorted(p[x] for x in t)) for t in ref} == ss for p in itertools.permutations(V)):
            iso += 1
    log(dict(stage='C0', labelled_complexes=len(sets), manifolds=len(good), genera=sorted({g for _, g in good}),
             isomorphic_to_torus7=iso, PASS=bool(good) and iso == len(good)))


# ---------------------------------------------------------------- Table 3 -> planes in torus7 numbering
def table3_planes():
    tri_of = {v: tuple(sorted(f for f in range(7) if v in T3F[f])) for v in T3V}
    assert all(len(t) == 3 for t in tri_of.values())
    target = set(S7.tris)
    sig = None
    for p in itertools.permutations(range(7)):
        if {tuple(sorted(p[x] for x in t)) for t in tri_of.values()} == target:
            sig = p; break
    assert sig is not None, 'Table 3 is not the dual of torus7'
    H = [None] * 7
    planar = True
    for f, vs in enumerate(T3F):
        P = [tuple(Q(x) for x in T3V[v]) for v in vs]
        nrm = None
        for a, b in itertools.combinations(range(1, 6), 2):
            nrm = cross(sub(P[a], P[0]), sub(P[b], P[0]))
            if any(nrm): break
        d = -dot(nrm, P[0])
        planar &= all(dot(nrm, p) + d == 0 for p in P)
        H[sig[f]] = (nrm[0], nrm[1], nrm[2], d)
        # the polygon order of Table 3 must equal the cyclic order given by the map
        k = sig[f]; L = S7.link[k]
        seq = [tuple(sorted(sig[x] for x in tri_of[v])) for v in vs]
        mseq = [tuple(sorted((k, L[s], L[(s + 1) % 6]))) for s in range(6)]
        ok = any(seq == mseq[r:] + mseq[:r] or seq[::-1] == mseq[r:] + mseq[:r] for r in range(6))
        assert ok, ('polygon order mismatch', f)
    return H, planar, sig


def chirotope_exact(H):
    rows = [tuple(Q(x) for x in h) for h in H] + [(Q(0), Q(0), Q(0), Q(1))]
    chi = {}
    for B in itertools.combinations(range(8), 4):
        d = check.det4([rows[b] for b in B])
        if d == 0:
            raise ValueError(('non-uniform', B))
        chi[B] = 1 if d > 0 else -1
    return chi


def corners_exact(H):
    """exact reflex bit per corner (k, frozenset(c, cp)) and convexity per face."""
    R, conv = {}, {}
    for k in range(7):
        L = S7.link[k]
        P = [check.point(H, k, L[s], L[(s + 1) % 6]) for s in range(6)]
        a = tuple(Q(x) for x in H[k][:3])
        area = sum(dot(cross(P[s], P[(s + 1) % 6]), a) for s in range(6))
        assert area != 0
        for s in range(6):
            t = dot(cross(sub(P[s], P[s - 1]), sub(P[(s + 1) % 6], P[s])), a)
            assert t != 0, 'straight corner'
            R[(k, frozenset((L[s], L[(s + 1) % 6])))] = (t > 0) != (area > 0)
        conv[k] = not any(R[(k, frozenset((L[s], L[(s + 1) % 6])))] for s in range(6))
    return R, conv


def enc1_eval(chi):
    bd = omsat.Builder(8, 'eval', chi)
    viol = omsat.dual_constraints(bd, S7, INF, ('orient',))
    return {k: len(v) for k, v in viol.items()}, dict(bd.Rloc), omsat.gp3_eval(chi, list(range(8)))


def enc2_reflex_lits(E):
    """corner (i, {j,k}) -> literal 'face i reflex at vertex ijk' (taken along L_ij)."""
    nv = E.nvar
    lits = {}
    for tri in S7.tris:
        for i in tri:
            j, k = [x for x in tri if x != i]
            lits[(i, frozenset((j, k)))] = E.reflex(i, j, k)
    assert E.nvar == nv, 'reflex literals were not all built already (cache miss)'
    return lits


def enc2_eval(chi):
    E = enc2_core.Enc(S7)
    chivals = np.array([[chi[B] > 0] for B in E.bases], dtype=bool)
    V = E.evaluate(chivals)
    fams = ('GP', 'SIMP', 'REFL', 'PIERCE', 'LOC')
    idx, mat = E.violated(V, fams)
    lits = enc2_reflex_lits(E)
    R = {c: bool(V[abs(l), 0] ^ (l < 0)) for c, l in lits.items()}
    return int(mat.sum()), R, E.defs_ok(V)


def perturb_uniform(H, Rref):
    """First seed whose exact perturbation is uniform, embedded and keeps all 42 reflex bits."""
    for seed in range(1000):
        rnd = random.Random(seed)
        Hp = [tuple(Q(x) + Q(rnd.randint(-5, 5), 10 ** 4) for x in h) for h in H]
        try:
            chi = chirotope_exact(Hp)
        except ValueError:
            continue
        okp, _ = check.check(S7, Hp)
        if okp and corners_exact(Hp)[0] == Rref:
            return seed, Hp, chi
    raise RuntimeError('no perturbation found')


def c1():
    H, planar, sig = table3_planes()
    ok, rep = check.check(S7, H)
    Rx, conv = corners_exact(H)
    seed, Hp, chi = perturb_uniform(H, Rx)
    log(dict(stage='C1-perturb', seed=seed, planes=[[str(x) for x in h] for h in Hp]))
    v1, R1, gpbad = enc1_eval(chi)
    v2, R2, defs2 = enc2_eval(chi)
    m1 = sum(R1[c] != Rx[c] for c in Rx)
    m2 = sum(R2[c] != Rx[c] for c in Rx)
    c = sum(conv.values())
    log(dict(stage='C1', table3_to_torus7=list(sig), planar=planar, checkerA_embedded=ok,
             checkerA_report={k: len(v) for k, v in rep.items()}, gp_violations=gpbad,
             enc1_violations=v1, enc1_reflex_mismatch=m1, enc2_violations=v2, enc2_defs_ok=defs2,
             enc2_reflex_mismatch=m2, corners=len(Rx), reflex_corners=sum(Rx.values()),
             convex_faces=c, convex_face_ids=[k for k in conv if conv[k]],
             PASS=planar and ok and gpbad == 0 and all(x == 0 for x in v1.values()) and m1 == 0 and v2 == 0 and m2 == 0
             and len(Rx) == 42))
    json.dump({str(B): s for B, s in chi.items()}, open(os.path.join(WORK, 'chi_table3.json'), 'w'))
    return chi, c


# ---------------------------------------------------------------- CNF writers
def card_clauses(C, m):
    """at least m of the literals C true: every (len(C)-m+1)-subset contains a true one."""
    if m <= 0:
        return []
    return [list(T) for T in itertools.combinations(C, len(C) - m + 1)]


def cnf_enc1(m, fix_chi=None, reorient=True):
    bd = omsat.Builder(8, 'cnf')
    omsat.gp3(bd, list(range(8)))
    omsat.dual_constraints(bd, S7, INF, ('orient',))
    if fix_chi is None and reorient:
        omsat.reorientation_fix(bd, list(range(8)))
    if fix_chi is not None:
        for B, s in fix_chi.items():
            bd.clauses.append([bd.var[B] if s > 0 else -bd.var[B]])
    C = []
    for k in range(7):
        ck = bd.new(); C.append(ck)
        corners = [key for key in bd.Rloc if key[0] == k]
        assert len(corners) == 6
        for key in corners:
            bd.clauses.append([-ck, -bd.Rloc[key]])
    bd.clauses += card_clauses(C, m)
    return bd.nv, bd.clauses, C


def cnf_enc2(m, fix_chi=None):
    E = enc2_core.Enc(S7)
    lits = enc2_reflex_lits(E)
    cls = E.cnf()
    nv = E.nvar
    if fix_chi is not None:
        for B, s in fix_chi.items():
            cls.append([E.bidx[B] if s > 0 else -E.bidx[B]])
    C = []
    for i in range(7):
        nv += 1; C.append(nv)
        corners = [c for c in lits if c[0] == i]
        assert len(corners) == 6
        for c in corners:
            cls.append([-nv, -lits[c]])
    cls += card_clauses(C, m)
    return nv, cls, C


def write(label, nv, cls):
    fn = os.path.join(WORK, 'cnf', label + '.cnf')
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn, 'w', newline='\n') as f:
        f.write(f'p cnf {nv} {len(cls)}\n')
        f.write(''.join(' '.join(map(str, c)) + ' 0\n' for c in cls))
    return fn


def random_bad_chi(seed=7):
    rnd = random.Random(seed)
    while True:
        H = [tuple(rnd.randint(-20, 20) for _ in range(4)) for _ in range(7)]
        try:
            chi = chirotope_exact(H)
            ok, _ = check.check(S7, H)
        except (ValueError, ZeroDivisionError):
            continue
        if not ok:
            return H, chi


def build_all():
    chi3 = {tuple(json.loads(k.replace('(', '[').replace(')', ']'))): v
            for k, v in json.load(open(os.path.join(WORK, 'chi_table3.json'))).items()}
    c = [json.loads(l) for l in open(OUT) if '"C1"' in l][-1]['convex_faces']
    Hbad, chibad = random_bad_chi()
    jobs = {
        'P1_enc1_m2': cnf_enc1(2),
        'P1b_enc2_m2': cnf_enc2(2),
        'P2_enc1_m1': cnf_enc1(1),
        'C4_enc1_base': cnf_enc1(0),
        'C4_enc2_base': cnf_enc2(0),
        'C3_enc1_badchi': cnf_enc1(0, fix_chi=chibad),
        'C3_enc1_m8': cnf_enc1(8),
    }
    if c < 7:
        jobs[f'C2_enc1_t3_m{c}'] = cnf_enc1(c, fix_chi=chi3)
        jobs[f'C2_enc1_t3_m{c + 1}'] = cnf_enc1(c + 1, fix_chi=chi3)
        jobs[f'C2_enc2_t3_m{c}'] = cnf_enc2(c, fix_chi=chi3)
        jobs[f'C2_enc2_t3_m{c + 1}'] = cnf_enc2(c + 1, fix_chi=chi3)
    meta = {}
    for lab, (nv, cls, C) in jobs.items():
        fn = write(lab, nv, cls)
        meta[lab] = dict(cnf=os.path.relpath(fn, WORK), vars=nv, clauses=len(cls), convex_vars=C,
                         sha256=hashlib.sha256(open(fn, 'rb').read()).hexdigest())
    json.dump(meta, open(os.path.join(WORK, 'cnf_meta.json'), 'w'), indent=1)
    log(dict(stage='CNF', bad_planes=[list(h) for h in Hbad], jobs={k: (v['vars'], v['clauses'], v['sha256'][:16])
                                                                    for k, v in meta.items()}))


def solve(label, limit=1800):
    meta = json.load(open(os.path.join(WORK, 'cnf_meta.json')))[label]
    fn = os.path.join(WORK, meta['cnf'])
    prf = os.path.join(WORK, 'proofs', label + '.drat')
    os.makedirs(os.path.dirname(prf), exist_ok=True)
    t0 = time.time()
    try:
        r = subprocess.run([CAD, fn, prf], capture_output=True, text=True, timeout=limit)
        res = {10: 'SAT', 20: 'UNSAT'}.get(r.returncode, f'rc{r.returncode}')
    except subprocess.TimeoutExpired:
        res, r = 'TIMEOUT', None
    t1 = time.time()
    rec = dict(stage='SOLVE', label=label, result=res, solve_s=round(t1 - t0, 2), cnf_sha256=meta['sha256'])
    if res == 'UNSAT':
        v = subprocess.run([DT, fn, prf, '-t', '40000'], capture_output=True, text=True)
        rec['drat_trim'] = 'VERIFIED' if 's VERIFIED' in v.stdout else 'CHECK-FAILED ' + v.stdout.strip()[-300:]
        rec['check_s'] = round(time.time() - t1, 2)
        rec['proof_bytes'] = os.path.getsize(prf)
        rec['proof_sha256'] = hashlib.sha256(open(prf, 'rb').read()).hexdigest()
    elif res == 'SAT':
        model = set()
        for line in r.stdout.splitlines():
            if line.startswith('v '):
                model |= {int(x) for x in line[2:].split() if x != '0'}
        rec['convex_true'] = [k for k, cv in enumerate(meta['convex_vars']) if cv in model]
        rec['chi_model'] = ''.join('+' if (q + 1) in model else '-' for q in range(70))
        if os.path.exists(prf):
            os.remove(prf)
    log(rec)
    return rec


def build_posthoc():
    """P1 without the reorientation fix (removes that normalisation from the argument)."""
    meta = json.load(open(os.path.join(WORK, 'cnf_meta.json')))
    for lab, (nv, cls, C) in {'PH_enc1_m2_nofix': cnf_enc1(2, reorient=False)}.items():
        fn = write(lab, nv, cls)
        meta[lab] = dict(cnf=os.path.relpath(fn, WORK), vars=nv, clauses=len(cls), convex_vars=C,
                         sha256=hashlib.sha256(open(fn, 'rb').read()).hexdigest())
    json.dump(meta, open(os.path.join(WORK, 'cnf_meta.json'), 'w'), indent=1)


def build_b():
    """Q instances: every face has a reflex corner (all 7 non-convex), plus their controls."""
    chi3 = {tuple(json.loads(k.replace('(', '[').replace(')', ']'))): v
            for k, v in json.load(open(os.path.join(WORK, 'chi_table3.json'))).items()}

    def e1(reorient=True, fix_chi=None):
        nv, cls, C = cnf_enc1(0, fix_chi=fix_chi, reorient=reorient)
        bd = omsat.Builder(8, 'cnf')               # rebuild to get the same Rloc literal numbering
        omsat.gp3(bd, list(range(8)))
        omsat.dual_constraints(bd, S7, INF, ('orient',))
        assert bd.nv == nv - 7 and C[0] == bd.nv + 1, 'literal numbering differs'
        for k in range(7):
            cls.append([bd.Rloc[key] for key in bd.Rloc if key[0] == k])
        return nv, cls, C

    def e2(fix_chi=None):
        nv, cls, C = cnf_enc2(0, fix_chi=fix_chi)
        lits = enc2_reflex_lits(enc2_core.Enc(S7))
        for i in range(7):
            cls.append([lits[c] for c in lits if c[0] == i])
        return nv, cls, C

    meta = json.load(open(os.path.join(WORK, 'cnf_meta.json')))
    jobs = {'Q1_enc1_nc7': e1(), 'Q2_enc2_nc7': e2(), 'Q3_enc1_nc7_nofix': e1(reorient=False),
            'D1_enc1_t3_nc7': e1(fix_chi=chi3), 'D1_enc2_t3_nc7': e2(fix_chi=chi3)}
    for lab, (nv, cls, C) in jobs.items():
        fn = write(lab, nv, cls)
        meta[lab] = dict(cnf=os.path.relpath(fn, WORK), vars=nv, clauses=len(cls), convex_vars=C,
                         sha256=hashlib.sha256(open(fn, 'rb').read()).hexdigest())
    json.dump(meta, open(os.path.join(WORK, 'cnf_meta.json'), 'w'), indent=1)
    log(dict(stage='CNF-B', jobs={k: meta[k]['sha256'][:16] for k in jobs}))


def build_d3():
    """Control: 'non-convex' required on six faces only must be SAT (the clauses are not contradictory by themselves)."""
    nv, cls, C = cnf_enc1(0)
    bd = omsat.Builder(8, 'cnf'); omsat.gp3(bd, list(range(8))); omsat.dual_constraints(bd, S7, INF, ('orient',))
    assert bd.nv == nv - 7
    for k in range(6):
        cls.append([bd.Rloc[key] for key in bd.Rloc if key[0] == k])
    meta = json.load(open(os.path.join(WORK, 'cnf_meta.json')))
    fn = write('D3_enc1_nc6', nv, cls)
    meta['D3_enc1_nc6'] = dict(cnf=os.path.relpath(fn, WORK), vars=nv, clauses=len(cls), convex_vars=C,
                               sha256=hashlib.sha256(open(fn, 'rb').read()).hexdigest())
    json.dump(meta, open(os.path.join(WORK, 'cnf_meta.json'), 'w'), indent=1)
