"""Phase 0: rescan the EXISTING near_*.npy snapshots (from the prior equal-weight GPU sweep) with the exact
checker and the genuine-piercing filter, to see whether any already beats the target (0 piercings, few crossed
faces) before spending new GPU time. Read-only re-use of existing data; writes only kp_ files."""
import glob, os, re, time
import numpy as np
import kp_common as kc

HERE = os.path.dirname(os.path.abspath(__file__))
results = []
files = sorted(glob.glob(os.path.join(HERE, 'near_*.npy')))
t0 = time.time()
for fp in files:
    fname = os.path.basename(fp)
    m = re.match(r'near_([A-D])_(\d+)', fname)
    if not m:
        m2 = re.match(r'near_seeded_(\d+)', fname)
        if not m2:
            continue
        num = int(m2.group(1))
    else:
        num = int(m.group(2))
    try:
        S = kc.surface_by_num(num)
    except KeyError:
        continue
    arr = np.load(fp)
    if arr.ndim == 2:
        arr = arr[None]
    nrows = min(3, arr.shape[0])
    for r in range(nrows):
        try:
            H = kc.to_rationals(arr[r])
            ok, rep, gX, sc = kc.full_report(S, H)
        except Exception as e:
            results.append((fname, r, num, None, str(e)))
            continue
        results.append((fname, r, num, sc, None))

results.sort(key=lambda x: (x[3] is None, x[3] if x[3] else (99, 99, 99)))
print(f'scanned {len(files)} files, {len(results)} rows, {time.time()-t0:.1f}s')
print(f"{'file':28s} {'row':3s} {'map':4s} {'score(pierce,faces,pts)':24s} err")
for fname, r, num, sc, err in results[:40]:
    errs = err or ''
    print(f'{fname:28s} {r:3d} {num:4d} {str(sc):24s} {errs}')

with open(os.path.join(HERE, 'kp_scan_existing_full.txt'), 'w') as f:
    for fname, r, num, sc, err in results:
        f.write(f'{fname}\t{r}\t{num}\t{sc}\t{err}\n')
