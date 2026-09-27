"""Reader's verification of this repository (run from the repository root).

  python reproduce.py hashes
      Regenerate every certified CNF from the encoder in code/ and compare its SHA-256 with results/certificates.csv.
      Needs numpy only. About 3 minutes.
  python reproduce.py prove MAP GROUP --cadical PATH --drat-trim PATH
      Regenerate one CNF (MAP as in the table, e.g. 54; GROUP e.g. T, D3, D2, S4, C3, C2, C2#1, C2#2), solve it with
      CaDiCaL writing a DRAT proof, and check the proof with drat-trim. The encoding is taken from the table.
      Example (7 s + 3 s):  python reproduce.py prove 54 T --cadical cadical --drat-trim drat-trim
  python reproduce.py controls
      Positive control: the 7-face (Szilassi) instance is satisfiable. Encoding against geometry: the clause evaluator
      agrees with the exact rational checker on random configurations. Needs numpy and python-sat.
"""
import argparse, csv, hashlib, os, subprocess, sys, time

ROOT = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(ROOT, 'code')
TABLE = os.path.join(ROOT, 'results', 'certificates.csv')


def rows():
    return list(csv.DictReader(open(TABLE)))


def extra_arg(enc):
    return {'plain': '-', 'orbfix': 'orbfix', 'orbfix+orient': 'orbfix,orient'}[enc]


def make_cnf(r):
    out = subprocess.run([sys.executable, 'write_cnf.py', r['map'], r['group'], extra_arg(r['encoding'])],
                         cwd=CODE, capture_output=True, text=True, check=True).stdout.strip()
    path = os.path.join(CODE, out)
    return path, hashlib.sha256(open(path, 'rb').read()).hexdigest()


def cmd_hashes(_):
    bad = 0
    for r in rows():
        path, h = make_cnf(r)
        ok = h == r['cnf_sha256']
        bad += not ok
        print(f"#{r['map']:>3} {r['group']:5} {r['encoding']:14} {'OK ' if ok else 'MISMATCH'} {h[:16]}", flush=True)
    print('ALL 36 CNFs REPRODUCED' if not bad else f'{bad} MISMATCHES')
    return bad


def cmd_prove(a):
    r = [x for x in rows() if x['map'] == a.map.lstrip('#') and x['group'] == a.group]
    if not r:
        sys.exit('no such case in the table')
    path, h = make_cnf(r[0])
    assert h == r[0]['cnf_sha256'], 'regenerated CNF differs from the certified one'
    proof = path[:-4] + '.drat'
    t0 = time.time()
    s = subprocess.run([a.cadical, '-q', path, proof], capture_output=True, text=True)
    print(f'solver exit code {s.returncode} (20 = UNSAT) after {time.time() - t0:.0f} s; proof {os.path.getsize(proof) / 1e6:.1f} MB')
    t0 = time.time()
    v = subprocess.run([a.drat_trim, path, proof], capture_output=True, text=True)
    line = [l for l in v.stdout.splitlines() if l.startswith('s ')]
    print(f"drat-trim: {line[0] if line else v.stdout[-300:]} after {time.time() - t0:.0f} s")
    return 0 if s.returncode == 20 and line and line[0].strip() == 's VERIFIED' else 1


def cmd_controls(_):
    sys.path.insert(0, CODE)
    os.chdir(CODE)
    import k12, run_sat
    from pysat.solvers import Solver
    S7 = k12.Surface(k12.torus7(), 't7')
    bd = run_sat.build('dual', S7)
    with Solver(name='cadical195', bootstrap_with=bd.clauses) as s:
        sat = s.solve()
    print('7-face (Szilassi) dual instance:', 'SAT (expected)' if sat else 'UNSAT (UNEXPECTED)')
    r = subprocess.run([sys.executable, 'validate_om.py', '0'], cwd=CODE, capture_output=True, text=True)
    print(r.stdout.strip())
    return 0 if sat and 'mismatches 0' in r.stdout else 1


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    sub.add_parser('hashes')
    q = sub.add_parser('prove'); q.add_argument('map'); q.add_argument('group')
    q.add_argument('--cadical', default='cadical'); q.add_argument('--drat-trim', dest='drat_trim', default='drat-trim')
    sub.add_parser('controls')
    a = p.parse_args()
    sys.exit({'hashes': cmd_hashes, 'prove': cmd_prove, 'controls': cmd_controls}[a.cmd](a))
