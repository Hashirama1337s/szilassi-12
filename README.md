# Szilassi polyhedra: two computer-assisted theorems

**Every ring of seven hexagons in which each face shares exactly one edge with every other has exactly one un-dented face, and the next possible such solid, with twelve faces, cannot have any symmetry: two computer-assisted proofs with proof files an independent checker reads.**

[![DOI 10.5281/zenodo.23003257](brand/doi.svg)](https://doi.org/10.5281/zenodo.23003257) ![Licence: text and data CC BY 4.0, code MIT](brand/licence.svg) [![Moki & Julio lab map](brand/lab.svg)](https://github.com/Hashirama1337s)

Moki & Julio · 2026 · version 1.1.1

**Read** [TECHNICAL.md](TECHNICAL.md) · **Part A** [convex7/README.md](convex7/README.md) · **Re-run** `python reproduce.py hashes`

## In plain words

Most solids you know have faces that touch only a few others. A cube's top face never touches its bottom. Only two
polyhedra are known in which **every** face shares exactly one edge with every other face: the tetrahedron, with 4 triangles,
and Szilassi's polyhedron (1977), a ring-shaped solid, like a lumpy bagel, made of 7 hexagons. (If some pairs of faces
may share two edges, two eight-faced examples are also known: Mizhaev 2020 and Röst and Vígh 2026.)

![Szilassi's polyhedron in three views: seven coloured hexagonal faces, a see-through hole, and the one un-dented face labelled](figures/szilassi-polyhedron.svg)

*Szilassi's polyhedron, drawn from the published exact coordinates: seven flat hexagons, one colour each. Every two of them share an edge (21 shared edges in all), and the solid is a ring with a hole you can see through. Six faces are dented, with one or two corners pointing inwards; face 7, the small face in the plane z = 4 (seen in the view from below), is the one un-dented face.*

Szilassi's hexagons are mostly dented: they have corners that point inwards. In 2009 Grünbaum and Szilassi conjectured
that every ring-shaped polyhedron built only from hexagons has at least six dented faces. **Part A** settles the
seven-hexagon case for rings in which no two faces are "overarching", that is, no two faces meet in more than one
place. This is a separate condition from each face sharing an edge with every other; for seven hexagons the two pick
out the same rings (see "Some words" below). The search covers every possible shape of that kind, not just Szilassi's model. Every one has exactly
one un-dented face: never two, never none.

**Part B** looks one step further. The next possible polyhedron in which every face shares exactly one edge with every other would have 12 faces.
Whether it exists has been open since the 1970s. We proved it cannot have any symmetry: no rotation, reflection or
other distance-preserving motion, apart from doing nothing, maps it onto itself. Whether a lopsided one exists is still open.

![](brand/icon-picture.svg) **The idea in one picture.** **Proving a sudoku has no solution.** You could try every filling, but a written chain of
deductions ("this square must be 3, so that one cannot be ...") ending in a contradiction is a proof anyone can
re-check line by line. Here the "squares" are yes/no facts about the shape (which side of each face's plane each corner
lies on), the "rules" say the faces may not cross each other, a SAT solver does the deducing, and it writes its chain
to a proof file that a separate small program re-checks.

**Where the picture breaks.** A sudoku's rules are exact. Here the rules are necessary conditions that every real solid
must meet, written by us as code, so the proof is only as good as that translation. That is why each translation was
tested against exact geometry, and for Part A written twice, independently. Extra solutions to the puzzle would be
harmless (they might not be real shapes); "no solution" is the answer that counts.

## At a glance

*The exact claims, for specialists. Terms are explained in "How it works" and "Words you will meet".*

- **Part A (v1.1, [`convex7/`](convex7/README.md)).** Let P be an acoptic polyhedral torus, in the sense of Grünbaum
  and Szilassi (2009, p. 21), with exactly seven faces, all hexagons, no two of them overarching. Then exactly one face
  of P is convex. This confirms their Conjecture 5.7 (p. 33), with equality, for seven-face hexangulations without
  overarching faces; those *with* overarching faces are not covered. Two independently written encodings (ENC1
  `code/omsat.py` with `orient`; ENC2 `code/enc2_core.py`) give UNSAT for "≥ 2 convex" and "no convex", every proof
  checked by drat-trim.
- **Part B (v1.0).** No polyhedron with 12 faces in which every two faces share exactly one edge has a non-trivial symmetry.
  Reduction (via the 59 neighbourly genus-6 triangulations and Aut(map)) to 36 SAT instances in an oriented-matroid
  (chirotope) encoding, all UNSAT with DRAT proofs verified by drat-trim. A second, independently derived encoding
  agrees on every case it completed; its half-turn runs did not finish at a 4-hour cap, so the half-turn cases rest on
  the main encoding's certificates alone.
- **Context.** The 59 maps are not regular, so Part B does not settle existence. Bokowski and H. (2025) conjecture
  that an embeddable regular map is also embeddable with a symmetry; under that heuristic (which does not
  literally cover these maps), the theorem is strong evidence that no 12-face Szilassi polyhedron exists.

Full technical record: [TECHNICAL.md](TECHNICAL.md).

## What we do not claim

- Seven-hexagon rings *with* overarching faces, and anything with eight or more faces (Part A).
- Whether a 12-face Szilassi polyhedron without symmetry exists: open (Part B). Long runs on several maps (up to 5 h
  plus 4 h) resolved none.
- Whether all 15 oriented-matroid classes found in Part A are realisable as actual solids.
- The SAT encodings use necessary conditions only; their correctness rests on the encoders, which are checked against
  exact geometry. For Part B's half-turn cases there is no second encoding.
- A referee reported, post hoc and without a certificate here, that the Part A upper bound may need fewer conditions,
  so a short hand proof may exist. Not checked.
- Correction 2026-11-17: "every face shares an edge with every other" means here "exactly one edge". Polyhedra in which
  some pairs of faces share two edges are not covered; two eight-faced examples of that kind are known (Mizhaev 2020;
  Röst and Vígh 2026).
- Correction 2026-11-17: "In plain words" called "no two faces overarching" the technical phrase for "each face shares
  an edge with every other". They are separate conditions; Part A is stated for the first, and for seven hexagons the
  two coincide (see "Some words").

## How it works

### The idea

```mermaid
flowchart LR
  A[Shape question:<br/>can such a solid exist?] --> B[Translate into yes/no facts<br/>+ rules every real solid obeys]
  B --> C[SAT solver:<br/>no way to satisfy all rules]
  C --> D[Proof file of the deduction]
  D --> E[drat-trim, a separate checker,<br/>re-checks every step: VERIFIED]
```

Some words. A face is **convex** (un-dented) if none of its corners point inwards. A corner that points inwards is
**reflex**. Two faces are **overarching** if they meet in more than one place; for seven hexagons, "no two overarching"
is the same as "every two share an edge".

### A worked example

**Part A in numbers.** Seven hexagons, each sharing an edge with each other one, give 7 × 6 / 2 = 21 edges. Euler's
formula for a ring (a torus) then forces 14 corners where faces meet, with exactly three faces at each. That pins
down one fixed pattern of faces and corners, Szilassi's map, so only the shape can vary. The computer asks two
questions about that map: "can at least two faces be convex?" and "can no face be convex?". Both come back UNSAT (no
solution), each with a proof that drat-trim verifies.

As a check, Szilassi's published coordinates, in exact fractions: all seven faces are flat, the solid is embedded (no
faces cross), 8 of its 42 face corners are reflex, and exactly one face, the one in the plane z = 4, is convex.

![The seven faces of Szilassi's polyhedron laid flat at the same scale, with inward corners circled: six faces have one or two, one face has none](figures/szilassi-faces-flat.svg)

*The seven faces laid flat, each in its own plane and all at the same scale. Open circles: corners that point inwards, 8 of the 42. Faces 1 to 6 each have one or two; face 7 has none. Part A proves that every seven-hexagon ring in which each face shares an edge with every other with no two faces "overarching" (a technical condition, see Part A; rings with overarching faces are not covered) has exactly one such un-dented face.*

**Part B in numbers.** A 12-face version would have twelve dented 11-sided faces, 44 corners and 66 edges, on a surface
with six holes (genus 6). There are exactly 59 possible patterns of how its faces could meet. A symmetry of the solid
would have to be a symmetry of its pattern, and that reduces the question to 25 cases (21 half-turns and 4 three-fold
rotations). With the larger groups as extra checks, that makes 36 SAT problems, all UNSAT. The fastest proof (the
12-element group T on pattern #54) takes 7 seconds; the half-turn cases take 3 to 135 minutes each.

![Dot chart of the 36 symmetric cases of the 12-face question, grouped by symmetry type, with the solver time for each impossibility proof on a log scale](figures/szilassi-symmetry-cases.svg)

*Part B: each dot is one possible symmetry acting on one of the 59 possible 12-face patterns, placed by how long the solver took to prove it impossible (log scale, 7 seconds to about 2 h 15 min). The 25 half-turn and one-third-turn cases settle the theorem; the larger groups are extra checks. All 36 proof files were re-checked by a separate checker: VERIFIED. A 12-face version without any symmetry is still open.*

## Reproduce and reuse

Questions this lets you answer:

- Do the proof instances really encode what they claim? (Rebuild every CNF and compare its SHA-256 with the table.)
- Is my own polyhedron embedded, with no faces crossing? (An exact checker reads any surface as JSON.)
- Can I attack the open case, lopsided 12-face polyhedra, myself? (The encoder writes those instances.)

**Rebuild and compare** (Python 3 with numpy; `controls` also needs python-sat). From the top folder:

```
git clone https://github.com/Hashirama1337s/szilassi-12
cd szilassi-12
python reproduce.py hashes        # ends "ALL 36 CNFs REPRODUCED"
python reproduce.py controls      # 7-face instance SAT; encoding vs geometry: 0 mismatches
cd convex7
python convex7.py checks          # seconds: ends "CHECKS PASS"
python convex7.py cnf             # ends "ALL 18 CNFs REPRODUCED"
```

To re-prove the instances themselves you also need CaDiCaL and drat-trim (on Windows, with the one-line binary-mode
patch in [`PATCHES.md`](PATCHES.md)). From the top folder:

```
python reproduce.py prove 54 T --cadical PATH --drat-trim PATH      # 10 seconds end to end
python reproduce.py prove 48 C2 --cadical PATH --drat-trim PATH     # about 6 minutes
cd convex7 && python convex7.py prove --cadical PATH --drat-trim PATH   # Part A, seconds (--full: 10-15 more min)
```

**Check your own polyhedron.** `code/checkerB.py` decides, in exact rational arithmetic, whether a polyhedral surface is
embedded. Input is JSON: `{"vertices": [["x", "y", "z"], ...], "faces": [[i0, i1, ...], ...]}`, with coordinates as
exact fractions written as strings.

```
python code/checkerB.py code/checkerB_tests/case08.json
# case08: EMBEDDED V=14 E=21 F=7 genus=1
```

On all 14 test polyhedra in `code/checkerB_tests/` it agrees with the separately written checker `code/check.py`.

**The open case.** `python code/write_cnf.py MAP full orient` writes the asymmetric instance for map number MAP (1 to
59) as a DIMACS CNF into `cnf/`, for anyone with a cluster or a new idea. A cube-based estimate puts these instances several
orders of magnitude beyond one PC.

## Words you will meet

- **Polyhedron** - a solid bounded by flat polygon faces.
- **Convex face** - a face with no corner pointing inwards; **reflex corner** - a corner that points inwards.
- **Torus / genus** - a ring shape; genus counts the holes (a bagel has 1; the 12-face case would have 6).
- **Overarching faces** - two faces that meet in more than one place.
- **Embedded** - built in real space with no faces crossing or piercing each other.
- **SAT solver** - a program that decides whether a large set of yes/no rules can all be satisfied at once.
- **UNSAT** - "no assignment satisfies all the rules".
- **DRAT proof / drat-trim** - a file recording the solver's deduction, and a small program that re-checks it.
- **CNF** - the standard text format for such a set of rules.
- **Symmetry** - a rotation, reflection or other distance-preserving motion that maps the solid onto itself.

## Open doors

What this work leaves open, and what we noticed along the way but did not follow up, for anyone who wants to take the next step.

- **A lopsided 12-face polyhedron.** Part B rules out every symmetric 12-face polyhedron in which each face shares exactly
  one edge with every other. Whether one with no symmetry at all exists is still unknown. *Status:* open according to
  Röst and Vígh 2026 (arXiv:2609.32998, Section 5: "the case of twelve faces, each pair sharing exactly one edge, remains
  open"). *Start here:* `code/write_cnf.py MAP full orient` writes the instance for any of the 59 face patterns;
  Section 7 of [TECHNICAL.md](TECHNICAL.md) lists what was tried, which encodings sped up the symmetric cases and which
  made the solver slower. Our long runs (up to 5 h plus 4 h) resolved none, and a cube-based estimate puts these instances
  several orders of magnitude beyond one PC. *Ceiling:* either a third known polyhedron in which every two faces share
  exactly one edge, or a proof that none with 12 faces exists.
- **The rest of the hexagon-ring conjecture.** Grünbaum and Szilassi conjectured that every ring-shaped polyhedron
  built only from hexagons has at least six dented faces. Part A settles seven-hexagon rings without overarching faces.
  Seven-hexagon rings *with* overarching faces, and rings of eight or more hexagons, are not covered. *Status:* open
  according to Grünbaum and Szilassi (2009), Conjecture 5.7, p. 33. *Start here:* [convex7/README.md](convex7/README.md),
  `convex7/convex7.py` and the two independent encodings (`code/omsat.py`, `code/enc2_core.py`). *Ceiling:* the
  conjecture settled in more cases, or a counterexample.
- **A short hand proof of "at most one un-dented face".** A referee reported, after the fact, that the "at most one
  convex face" half of Part A stays unsolvable when the edge-piercing and one-dented-corner-per-vertex rules are removed.
  If so, it needs only that faces do not cross themselves and the sign rules of the seven face planes, and a short
  pen-and-paper proof may exist. *Status:* UNTESTED observation from this work (reported without a certificate here).
  *Start here:* `convex7/convex7.py cnf` rebuilds the instances; remove the clause families named above and re-run
  with a DRAT checker. *Ceiling:* a proof a reader can check by hand for half of Part A.
- **Which of the 15 sign patterns are real solids?** Without any convexity rule, the main encoding of the seven-hexagon
  ring has exactly 1,260 solutions, which fall into 15 classes up to relabelling and mirror image. They differ in how
  the plane at infinity sits among the seven face planes. Whether every one of the 15 is realised by an actual solid
  is not decided. *Status:* UNTESTED observation from this work. *Start here:* `convex7/om_classes.py` (the count) and
  the exact embeddedness checker `code/checkerB.py`. *Ceiling:* a complete list of the different ways a seven-hexagon
  ring of this kind can sit in space, at the level of these sign patterns.

## How to cite

```
Moki & Julio (2026). Szilassi polyhedra: every seven-hexagon torus without overarching faces has exactly one convex face, and no symmetric 12-face Szilassi polyhedron exists (computer-assisted proofs with DRAT certificates). Version 1.1.1. Zenodo. https://doi.org/10.5281/zenodo.23003257
```

The DOI always resolves to the newest archived version. Version 1.0 (Part B only) is 10.5281/zenodo.23003258. Machine-readable: [CITATION.cff](CITATION.cff).
Text and data are under CC BY 4.0 (`LICENSE`); code is under the MIT licence (`LICENSE-CODE`).

Release note - v1.1.1: adds Open doors and figures; corrections dated 2026-11-17; no change to the certified results.

## Related projects

- [moki-julio-circle-packing](https://github.com/Hashirama1337s/moki-julio-circle-packing): record packings of equal circles, spheres and hyperspheres, with exact certificates
- All projects: [the Moki & Julio lab map](https://github.com/Hashirama1337s)

---

*Prepared with AI assistance; what was checked, and by what, is stated above.*
