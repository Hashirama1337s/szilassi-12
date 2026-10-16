# Szilassi polyhedra: two computer-assisted theorems

**Computer-assisted theorems with checkable proof certificates. Moki & Julio (AI-assisted), 2026.**

## In plain words

Szilassi's polyhedron is a ring-shaped solid made of seven flat hexagons, and each hexagon shares an edge with every other one.
This repository proves two facts about polyhedra of this kind:

- **Part A (v1.1, folder `convex7/`):** every seven-hexagon ring in which every two faces share an edge has exactly
  one un-dented (convex) face, never two and never none. This settles a 2009 conjecture of Grünbaum and Szilassi for
  seven-face rings of that kind.
- **Part B (v1.0):** the next possible "every face touches every other" polyhedron would have twelve faces. It cannot
  have any symmetry.

Each "impossible" answer comes with a proof file checked by an independent checker program.

## Part A: every seven-hexagon torus without overarching faces has exactly one convex face

> **Theorem.** Let P be an acoptic polyhedral torus, in the sense of Grünbaum and Szilassi (*Contrib. Discrete Math.*
> 4(1) (2009) 21–39, p. 21), with exactly seven faces, all hexagons, no two of them overarching. Then exactly one
> face of P is convex.

**What it settles.** Grünbaum and Szilassi conjectured (Conjecture 5.7, p. 33) that every hexangulation of the torus
has at least six non-convex faces. The theorem confirms this, with equality, for all seven-face hexangulations without
overarching faces. Seven-face hexangulations *with* overarching faces are not covered.

**How it was proved.**

- Two independently written SAT encodings of embedded realisations, plus convex-face constraints, give UNSAT for "at
  least two convex faces" and for "no convex face".
- Each UNSAT answer has a DRAT proof verified by drat-trim.
- Szilassi's published coordinates pass every check and have exactly one convex face.

Details, controls and a one-command reproduction are in [`convex7/README.md`](convex7/README.md).

## Part B: no symmetric 12-face Szilassi polyhedron

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23003257.svg)](https://doi.org/10.5281/zenodo.23003257) DOI 10.5281/zenodo.23003257 (all versions; v1.0: 10.5281/zenodo.23003258)

Only two polyhedra are known in which every two faces share an edge: the tetrahedron (4 faces) and Szilassi's
polyhedron (7 faces, genus 1, 1977). The next possible case has 12 faces. It would be a polyhedron of genus 6 whose
faces are twelve non-convex 11-gons, with 44 vertices and 66 edges. Whether it exists has been open since the 1970s
(Wikipedia, "Szilassi polyhedron"; Bokowski and H., *Symmetry* 17 (2025) 622, sec. 8.3).

> **Theorem.** No polyhedron with 12 faces in which every two faces share an edge has a non-trivial symmetry.
> No rotation, rotoreflection or other isometry other than the identity maps such a polyhedron onto itself.

The proof is a reduction to 36 SAT instances, all unsatisfiable. Each instance has a DRAT proof checked by
[drat-trim](https://github.com/marijnheule/drat-trim), and an independently derived second encoding agrees on every case
it completed. The question for asymmetric polyhedra remains open. Section 7 says what was tried and how hard it looks.

Bokowski and H. (2025, sec. 4) report that every embedding they found had a non-trivial symmetry. They conjecture
that an embeddable regular map is also embeddable with one. The 59 maps below are not regular, so the theorem does not
settle existence. Under that heuristic, it is strong evidence that no 12-face Szilassi polyhedron exists.

## 1. Setting

Combinatorially, such a polyhedron is the dual of a triangulation of the orientable genus-6 surface with 12 vertices in
which every two vertices are joined by an edge. There are exactly 59 such triangulations (Altshuler, Bokowski and
Schuchert, *J. Combin. Theory A* 75 (1996)). We use F. H. Lutz's list `manifold_lex_d2_n12_o1_g6` and his numbering
#1 to #59 (`code/k12_g6_lutz.txt`, sha256 `55912614618733476d07094ca4ed638df84ab58b42e8dba9eb81730da1a589b7`).

Schewe (*Discrete Comput. Geom.* 43 (2010); arXiv:0801.2582) showed with oriented matroids and SAT solvers that the
*vertex* version (12 vertices, all pairs joined) cannot be realised. To our knowledge the face version had not been
attacked this way. Its only published computer attempt is the numerical search of Bokowski and H. (2025), which found
no embedding.

## 2. Reduction to 25 cases

- 20 of the 59 maps have non-trivial automorphisms and 39 have none (`code/autos.py`, `code/automorphisms.json`).
- A symmetry of a realisation permutes the faces. Only the identity fixes all 12 face planes, so the symmetry group
  embeds in Aut(map). Every non-trivial automorphism of every map fixes no face, so the group acts freely on the faces.
- The element orders that occur are 2, 3 and 4. Elements of order 2 preserve orientation; elements of order 4 reverse
  it and square to order 2. So every non-trivial symmetry group contains a half-turn or a 3-fold rotation. A finite
  isometry group fixes a point, which we place at the origin.
- Up to conjugation in Aut(map) this leaves 21 half-turn cases (19 maps; #59 has three classes of involutions) and 4
  three-fold cases (#3, #45, #48, #54). All half-turns are conjugate in SO(3). The two directions of a 3-fold rotation
  give mirror images of the same unoriented map.
- The larger groups (T on #54; D3 on #3 and #48; D2 on #54 and #59; the 4-fold rotoreflection S4 on #5, #6, #23, #41,
  #57 and #58) follow from these 25 cases. They were proved separately as extra checks.

## 3. The encoding

A configuration is 12 planes h_i = (a_i, d_i) (plane a.x + d = 0) together with the plane at infinity h_INF = (0,0,0,1).
Its rank-4 chirotope chi assigns each 4-element subset of the 13 the sign of its determinant (715 signs). The polyhedron's
vertices are the points v(i,j,k) where three planes meet. With chi(i,j,k,INF) = sign det(a_i,a_j,a_k):

- **Order along a line.** On L_ij (direction a_i x a_j), v(i,j,l) comes after v(i,j,k) iff
  chi(ijlk) chi(ijl INF) chi(ijk INF) = +.
- **Point on an edge.** v(i,j,k) lies strictly inside the edge ij, whose ends are v(i,j,a) and v(i,j,b), iff
  chi(ijka) chi(ija INF) chi(ijkb) chi(ijb INF) = -.

The polyhedron is embedded iff **(S)** no face crosses itself and **(X)** no edge meets a face it does not belong to.
Every two faces are adjacent, so every contact between faces reduces to these two cases. "v(i,j,k) is inside face k" is
the parity of boundary crossings along a ray on L_ki. The 3-term Grassmann-Pluecker relations, which characterise
uniform chirotopes, complete the CNF (`code/omsat.py`). Only necessity is used: every embedded realisation gives a
satisfying assignment. A symmetric instance adds chi(gB) = det(g) chi(B).

**Orientation encoding** (`orient`, used for the #23 and #48 half-turn certificates). Each face gets a winding bit
W_k. A corner of face k is reflex exactly when its local turn disagrees with W_k. Two constraints are added:

1. the reflex bit computed by parity equals turn XOR W_k at every corner;
2. at most one reflex corner meets at any polyhedron vertex. Two reflex corners at a degree-3 vertex would make those
   faces overlap along their common line.

Both hold in every embedded realisation. This made the half-turn instances about 10x faster.

**Steps that lose no generality.**

1. *General position.* The realisation space is open. In all 45 (map, group) cases, a random symmetric configuration
   has all 715 determinants non-zero (`code/uniform_check.py`), so any symmetric realisation can be perturbed
   symmetrically into general position.
2. *Reorientation.* Replacing h by -h on a whole face orbit, or on h_INF, changes no plane and no constraint. One sign
   is fixed per independent flip class (`orbfix`, 12x faster).

## 4. Verification

- **Controls.**
  - The 7-face instance (Szilassi) is satisfiable.
  - The same GPU search re-finds Szilassi polyhedra from random starts.
  - The 7-vertex Csaszar instance (vertex version) is satisfiable.
  - Schewe's theorem is reproduced: all 59 vertex versions are UNSAT, about 4 minutes each (`results/logs/sat_schewe.log`).
- **Encoding against geometry.** On 45 configurations (random and near-miss), the clause evaluator agrees test-by-test
  with an exact rational geometric checker, with 0 mismatches:
  - 17,715 face-crossing tests, 1,307 of them violated;
  - 10,214 edge-face tests, 431 of them violated.

  The orientation constraints give 0 violations on real Szilassi polyhedra (`code/validate_*.py`).
- **Independent second encoding** (`code/enc2_*.py`). It was written from the geometric definition without access to
  the first, using a different primitive: which side of a plane a triple point lies on.
  - It matches its own exact geometric test on 104 configurations (102,417 conditions, 0 mismatches).
  - It finds the 7-face case satisfiable.
  - It proves #54 T, #3 D3, #48 D3, #59 D2, #5 S4 and #45 C3 unsatisfiable.
  - Its half-turn runs (#55, #59) were still unresolved at a 4-hour cap on a shared machine. So the half-turn cases
    rest on the main encoding's certificates alone.
- **Two exact embeddedness checkers**, written separately:
  - `code/check.py`: edge-face tests in rational arithmetic;
  - `code/checkerB.py`: ear clipping plus exact triangle-triangle tests.

  They agree on all 14 test polyhedra in `code/checkerB_tests/`.
- **Proof certificates.** Every one of the 36 instances was solved by CaDiCaL writing a DRAT proof, which drat-trim
  verified. In total the proofs came to about 8.4 GB; each is regenerated on demand rather than stored. The table
  `results/certificates.csv` gives, for each case, the encoding, solver, times, proof size, the verdict (all
  VERIFIED) and the full SHA-256 of the CNF. Raw logs are in `results/logs/certified_*.log`.

| group | maps | solve time | proof |
|---|---|---|---|
| T (12 elements) | #54 | 7 s | 0.1 MB |
| D3 | #3, #48 | 26-31 s | 35 MB |
| D2 | #54, #59 | 134-138 s | 45-46 MB |
| S4 | #5, #6, #23, #41, #57, #58 | 186-281 s | 46-68 MB |
| C3 | #3, #45, #48, #54 | 170-255 s | 47-64 MB |
| C2 (half-turn) | 19 maps, 21 classes | 3-135 min | 42-691 MB |

## 5. Reproduce

Requirements: Python 3 with numpy; python-sat for `controls`; CaDiCaL and drat-trim for `prove`.

```
python reproduce.py hashes                       # regenerate all 36 CNFs, compare SHA-256 with the table
python reproduce.py prove 54 T --cadical PATH --drat-trim PATH      # 10 seconds end to end
python reproduce.py prove 48 C2 --cadical PATH --drat-trim PATH     # about 6 minutes
python reproduce.py controls                     # 7-face SAT + encoding vs geometry
```

On Windows, both tools need a one-line fix to open proof files in binary mode (see `PATCHES.md`). Other scripts:
`code/sym_sat.py MAP GROUP [orbfix[,orient]]` solves one symmetric instance, and `code/full_sat.py MAP [orient]` solves
a full (asymmetric) instance.

## 6. Numerical side

A batched GPU search (`code/search.py`) minimised the violation count used by Bokowski and H. (edge-face hits plus face
self-crossings). It found no realisation.

- Symmetric searches got as low as 4 violations (S4 on #5 and #58); those classes are now proved impossible.
- Unrestricted searches over all 59 maps plateau at 21-24.
- A dedicated search for the self-intersecting ("Kepler-Poinsot") version with the fewest crossings found map #5 with
  no edge piercing a face and four faces crossing themselves once each. This was verified exactly
  (`results/near_miss/`). It does not beat the published model with two crossed polygons.

![near miss](results/near_miss/best_5.png)

## 7. Open: asymmetric realisations

The asymmetric case covers the 39 asymmetric maps, and asymmetric realisations of the other 20.

- **Size.** The full instances have 715 free signs, about twice as many as the half-turn ones.
- **Long runs.** None was resolved: #1 and #54 ran 5 h with the base encoding and then 4 h with the orientation
  encoding; #5 and #59 ran 4 h with the orientation encoding; the other 55 maps ran 10 min each.
- **Encodings that did not help.** Four redundant encodings made the solver slower: explicit order transitivity,
  Jordan consistency between two rays, the 3-D even-crossing rule on lines offset from the edges, and propagating
  inside-status between consecutive points.
- **Encodings that did help.** Reorientation fixing (12x) and the orientation encoding (2.7x on 3-fold cases, 10x on
  half-turn cases, growing with size).
- **No local lemma.** The proof cores are global: the #45 C3 proof uses 503 of the 528 face-simplicity constraints and
  405 of the 528 edge-face constraints.

A cube-based estimate puts the asymmetric instances several orders of magnitude beyond one PC. The encoder generates
them (`code/write_cnf.py MAP full orient`), for anyone with a cluster or a new idea.

## Credits and licence

- Data: the 59 triangulations are from Altshuler, Bokowski and Schuchert (1996), as listed on F. H. Lutz's Manifold Page.
- Part A uses the coordinates of Szilassi's polyhedron from Grünbaum and Szilassi (2009), Appendix Table 3
  (doi:10.55016/ojs/cdm.v4i1.61986).
- v1.1 adds one line to `code/omsat.py`: it exposes the per-corner reflex literals for `convex7/`. The line changes no
  CNF, and `python reproduce.py hashes` still reproduces all 36 Part B certificates.
- Tools: CaDiCaL (A. Biere et al.), drat-trim (M. Heule et al.), PySAT.
- Text and data are under CC BY 4.0 (`LICENSE`); code is under the MIT licence (`LICENSE-CODE`). Please cite via
  `CITATION.cff`.
- Two staged copies differ from the working files in small ways: `code/symgroup.py` imports torch lazily (the numerical
  search needs it, CNF generation does not), and `code/make_tests.py` writes its answer key to
  `checkerB_tests/expected_checkerA.json`.

## Other work by Moki & Julio

- [moki-julio-wifi-ris-cell](https://github.com/Hashirama1337s/moki-julio-wifi-ris-cell): One 1-bit reconfigurable-intelligent-surface cell covering the 2.4 GHz and 5-7 GHz Wi-Fi bands with one switch state (simulation study, Palace FEM) ([doi:10.5281/zenodo.23165406](https://doi.org/10.5281/zenodo.23165406))
- [rikitake-chaos](https://github.com/Hashirama1337s/rikitake-chaos): Computer-assisted proof that Rikitake's two-disc dynamo (1958) is chaotic ([doi:10.5281/zenodo.23041182](https://doi.org/10.5281/zenodo.23041182))
- [moki-julio-circle-packing](https://github.com/Hashirama1337s/moki-julio-circle-packing): 6,071 new best-known packings of equal circles, spheres and hyperspheres, each verified by two independent exact checkers ([doi:10.5281/zenodo.22981305](https://doi.org/10.5281/zenodo.22981305))
- [kinect-v1-depth-levels](https://github.com/Hashirama1337s/kinect-v1-depth-levels): The Xbox 360 Kinect measured as an instrument: the depth lattice, three noise populations, drift and a capability map

All projects: [github.com/Hashirama1337s](https://github.com/Hashirama1337s)
