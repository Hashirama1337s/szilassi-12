"""Reader's verification of the 7-face convexity theorem (run from this folder, convex7/).

  python convex7.py checks
      C0: the 7-vertex torus triangulation is unique (all 120 labelled ones are isomorphic).
      C1: Szilassi's coordinates (Grunbaum-Szilassi 2009, Table 3): planar, embedded (exact checker), exactly one convex
          face; after an exact perturbation into general position both encodings show no violation and their reflex bits
          match exact geometry at all 42 corners.
  python convex7.py cnf
      Writes every instance to out/cnf/ and compares its SHA-256 with certificates.csv.
  python convex7.py prove [LABEL ...] [--full]
      Solves the instances with CaDiCaL (writing a DRAT proof) and checks each UNSAT proof with drat-trim; compares every
      result with certificates.csv. Without labels: every instance except the four slow ones (add --full for those,
      about 10-15 minutes more, mostly drat-trim). Executables: --cadical PATH --drat-trim PATH (or env CADICAL / DRAT_TRIM).
Needs Python 3 with numpy."""
import argparse, csv, json, os, sys, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = os.path.join(HERE, 'certificates.csv')
SLOW = {'P1b_enc2_m2', 'PH_enc1_m2_nofix', 'Q2_enc2_nc7', 'Q3_enc1_nc7_nofix'}


def table():
    return {r['label']: r for r in csv.DictReader(open(TABLE, encoding='utf-8'))}


def cmd_checks(core):
    core.c0()
    core.c1()
    recs = [json.loads(l) for l in open(core.OUT, encoding='utf-8')]
    ok = [r for r in recs if r['stage'] == 'C0'][-1]['PASS'] and [r for r in recs if r['stage'] == 'C1'][-1]['PASS']
    print('CHECKS', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


def cmd_cnf(core):
    if not os.path.exists(os.path.join(core.WORK, 'chi_table3.json')):
        core.c1()
    core.build_all()
    core.build_posthoc()
    core.build_b()
    core.build_d3()
    meta = json.load(open(os.path.join(core.WORK, 'cnf_meta.json')))
    bad = 0
    for lab, r in table().items():
        got = meta[lab]['sha256']
        bad += got != r['cnf_sha256']
        print(f"{lab:22} {'OK ' if got == r['cnf_sha256'] else 'MISMATCH'} {got[:16]}")
    print(f'ALL {len(table())} CNFs REPRODUCED' if not bad else f'{bad} MISMATCHES')
    return 1 if bad else 0


def cmd_prove(core, labels, full):
    if not os.path.exists(os.path.join(core.WORK, 'cnf_meta.json')):
        if cmd_cnf(core):
            return 1
    T = table()
    labels = labels or [l for l in T if full or l not in SLOW]
    bad = 0
    for lab in labels:
        rec = core.solve(lab)
        want = T[lab]['result']
        ok = rec['result'] == want and (want != 'UNSAT' or rec.get('drat_trim') == 'VERIFIED')
        bad += not ok
        print(f"{lab:22} expected {want:5} got {rec['result']:5} {rec.get('drat_trim', '')} {'OK' if ok else 'FAIL'}")
    print('ALL RESULTS REPRODUCED' if not bad else f'{bad} FAILURES')
    return 1 if bad else 0


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=['checks', 'cnf', 'prove'])
    p.add_argument('labels', nargs='*')
    p.add_argument('--full', action='store_true')
    p.add_argument('--cadical'); p.add_argument('--drat-trim', dest='drat_trim')
    a = p.parse_args()
    if a.cadical:
        os.environ['CADICAL'] = a.cadical
    if a.drat_trim:
        os.environ['DRAT_TRIM'] = a.drat_trim
    sys.path.insert(0, HERE)
    import core
    sys.exit({'checks': lambda: cmd_checks(core), 'cnf': lambda: cmd_cnf(core),
              'prove': lambda: cmd_prove(core, a.labels, a.full)}[a.cmd]())
