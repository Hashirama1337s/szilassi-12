"""Certify UNSAT instances: CaDiCaL (standalone) writes a DRAT proof, drat-trim checks it.
usage: python certify.py TAG KEY:LABEL [KEY:LABEL ...]   -> certified_<TAG>.log ; proofs deleted after a VERIFIED check
(the CNF + solver + checker reproduce them); failures keep their proof."""
import sys, os, subprocess, time, hashlib
tag, jobs = sys.argv[1], sys.argv[2:]
log = open(f'certified_{tag}.log', 'a', encoding='utf-8')
CAD, DT = r'tools\cadical195\cadical.exe', r'tools\drat-trim\drat-trim.exe'   # CaDiCaL 1.9.5 (patched 'wb'); drat-trim patched 'rb'
os.makedirs('proofs', exist_ok=True)
for job in jobs:
    key, lab, *ex = job.split(':')
    cnf = subprocess.run([sys.executable, 'write_cnf.py', key, lab] + ex, capture_output=True, text=True).stdout.strip()
    sha = hashlib.sha256(open(cnf, 'rb').read()).hexdigest()[:16]
    prf = os.path.join('proofs', os.path.basename(cnf).replace('.cnf', '.drat'))
    t0 = time.time()
    r = subprocess.run([CAD, '-q', cnf, prf], capture_output=True, text=True)
    res = {10: 'SAT', 20: 'UNSAT'}.get(r.returncode, f'rc{r.returncode}')
    t1 = time.time()
    if res == 'UNSAT':
        v = subprocess.run([DT, cnf, prf, '-t', '40000'], capture_output=True, text=True)
        ok = 's VERIFIED' in v.stdout
        size = os.path.getsize(prf)
        verdict = 'VERIFIED' if ok else 'CHECK-FAILED: ' + v.stdout.strip()[-200:]
        if ok:
            os.remove(prf)
    else:
        verdict, size = '-', 0
    log.write(f'{time.ctime()} #{key} {lab}: {res} solve {t1 - t0:.0f}s proof {size / 1e6:.0f} MB drat-trim {time.time() - t1:.0f}s {verdict} cnf {os.path.basename(cnf)} sha256 {sha}\n')
    log.flush()
log.write('DONE\n')
