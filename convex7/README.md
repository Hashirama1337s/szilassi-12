# Every seven-hexagon torus without overarching faces has exactly one convex face

**A computer-assisted theorem with checkable proof certificates. Moki & Julio (AI-assisted), 2026.**

## In plain words

Szilassi's polyhedron is a ring-shaped solid made of seven flat hexagons, and each hexagon shares an edge with every other one. In
2009 Grünbaum and Szilassi conjectured that every ring-shaped polyhedron built only from hexagons has at least six
dented (non-convex) faces.

We settle the seven-hexagon case for rings in which every two hexagons share an edge, as in Szilassi's model. The
computer search covers every possible shape of that kind, not just Szilassi's model. In every such polyhedron
exactly one face is un-dented: never two, and never none. So for these rings six dented faces is both forced and the
minimum, as conjectured. Rings in which two faces meet in more than one place are not covered.

Each "impossible" answer comes with a proof file that an independent checker program verifies. Two separately written
encodings of the geometry give the same answers.

## 1. Statement

We use the definitions of Grünbaum and Szilassi [GS, p. 21]:

- a face is a planar simple polygon, with no collinear adjacent edges;
- adjacent faces are not coplanar;
- the polyhedron is *acoptic*: two faces meet only in common vertices or edges;
- two faces are *overarching* if their intersection has more than one component;
- the faces at each vertex form one circuit, so at least three faces meet there [GS, p. 21].

A face is *convex* if it is a convex polygon.

> **Theorem.** Let P be an acoptic polyhedral torus with exactly seven faces, all hexagons, no two of them
> overarching. Then exactly one face of P is convex.

**Corollary.** Conjecture 5.7 of [GS, p. 33] says "every hexangulation of the torus has at least six non-convex faces".
It holds, with equality, for every seven-face hexangulation without overarching faces.

**Scope.**

- The conjecture also covers hexangulations with overarching faces, which [GS, Remark 5.4] admits. Seven-face ones of
  that kind are **not** covered here.
- Nothing is claimed for eight or more faces.
- For seven hexagons, "no two faces overarching" is the same as "every two faces share an edge". If every pair
  shares an edge, then E = 21 forces exactly one shared edge per pair and three faces at each vertex, so two faces
  meet only along that edge.
- Szilassi's model [GS, Table 3] has exactly one convex face (checked exactly, §5), as the theorem requires.

**Prior work.**
- That a hexangulation cannot have only convex faces is called well known in [GS, Remark 5.6, p. 33].
- For seven faces, "at most four convex faces" (so at least three non-convex) already follows from published
  work. Every two faces of P share an edge, and no surface of convex polygons has five polygons that pairwise share
  sides ([AK, §3.1, Prop. 6], credited there to Thomassen 1997).
- Particular models with exactly one convex face are described in [EJ] and [MC].
- Bokowski and Schewe (2002) are reported to have shown, with oriented matroids, that Szilassi's torus is
  essentially unique: "axially symmetric" in [SV], unqualified in [BH, §8.2]. We could not obtain their paper.
  Uniqueness of the oriented matroid of the seven face planes does not by itself decide convexity, which also
  depends on where the plane at infinity lies (see §5, "Oriented matroids").
- We found no published statement of the theorem, or of any bound on the convex faces beyond the above.

## 2. Reduction to one map (by hand)

1. **Every two faces share an edge.**
   - Two faces that are not overarching meet in at most one edge.
   - A hexagon has six edges, so each face has six different neighbours: all the other faces.
2. **Every vertex lies on exactly three faces.**
   - Euler's formula for the torus with F = 7 and E = 21 gives V = 14.
   - Every vertex lies on at least three faces [GS, p. 21], and the 42 edge-ends spread over 14 vertices give exactly
     three per vertex.
3. **The dual is the Möbius torus.**
   - The dual is a triangulated torus on seven vertices in which every two vertices are joined.
   - It is a simplicial complex. If two vertices of P lay on the same three faces, the three edges between those faces
     would all join the same two points and coincide.
   - The seven-vertex triangulation of the torus is unique. `convex7.py checks` enumerates all 120 labelled ones and
     finds them all isomorphic to triangles {i, i+1, i+3} and {i, i+2, i+3} (mod 7).
   - So P realises one fixed map: the dual of that triangulation, Szilassi's map.
4. **General position (sketch).**
   - The three planes at each vertex have independent normals. Otherwise the three edges at that vertex would lie on
     one line and two of them would overlap.
   - So the 14 vertices are simple triple points, and they move continuously when the seven planes move slightly.
   - Under a small motion, each face stays a simple hexagon with the same reflex corners. No corner is straight
     (no collinear adjacent edges [GS, p. 21]), and non-adjacent edges of a face are a positive distance apart.
   - The polyhedron stays embedded. Features that do not meet (an edge and a face not containing it, two faces away
     from their shared edge) are disjoint compact sets at positive distance. The local picture at every edge and
     vertex (two or three planes in general position) is stable.
   - The configurations in which some 4×4 determinant of the planes (with the plane at infinity added) vanishes form
     a finite union of hypersurfaces, which a small motion avoids.

   So every realisation can be moved to one with a uniform chirotope, keeping it embedded and keeping its convex-face
   count. Strict and non-strict convexity coincide throughout.

## 3. The SAT encodings

The variables are the signs of the 70 determinants (a rank-4 chirotope on 8 elements). Two encodings, written
independently, each turn the necessary conditions for an embedded realisation into clauses:

- **ENC1:** `code/omsat.py` with its `orient` option. It is the main encoding of this repository (Part B, §3 of the
  top-level README).
  - It includes the 3-term Grassmann–Plücker relations, face simplicity and edge-face non-piercing.
  - Each face has a winding bit, and each corner a reflex bit that must agree with a parity test along the edge line.
  - At most one reflex corner meets at any vertex.
- **ENC2:** `code/enc2_core.py`, written separately, from the geometric definition. Its conditions are simplicity,
  piercing, a reflex predicate and its local consistency, plus Grassmann–Plücker.

**Convexity.** Face k gets an indicator C_k, with clauses C_k → (corner not reflex) for each of its six corners.

- "At least two convex faces" adds the 7 clauses OR(C_k : k ∈ T), one for each 6-element set T of faces.
- "No convex face" instead adds, for each face, the clause "some corner of this face is reflex".

Only necessity is used: every embedded general-position realisation gives a satisfying assignment, with its true signs
and true reflex bits. So UNSAT means no such polyhedron exists.

**What is trusted.**
- The DRAT proofs certify the CNF formulas.
- That the formulas encode necessary conditions rests on the two encoders. They were written independently and are
  checked against exact geometry (§5).
- The search space is uniform sign patterns satisfying the 3-term Grassmann–Plücker relations. This is a relaxation
  of realisability, which is harmless for UNSAT answers.

One normalisation is optional. Reversing a plane's orientation changes no polyhedron, so ENC1 can fix one sign per
reorientation class. Every theorem instance is also solved without it.

## 4. Results

All runs used CaDiCaL 1.9.5, which wrote a DRAT proof, and drat-trim, which checked it. The full table, with SHA-256
of every CNF and proof, is in `certificates.csv`.

| question | ENC1 + fix | ENC1, no fix | ENC2 |
|---|---|---|---|
| at least 2 convex faces | UNSAT 0.46 s, VERIFIED | UNSAT 48 s, VERIFIED | UNSAT 47 s, VERIFIED |
| no convex face | UNSAT 0.36 s, VERIFIED | UNSAT 24 s, VERIFIED | UNSAT 51 s, VERIFIED |
| at least 1 convex face | SAT | – | – |

The six theorem proofs are 0.57–32 MB.

## 5. Controls

Every control comes out as expected.

- **C0:** all 120 labelled seven-vertex torus triangulations are isomorphic (step 3 above).
- **C1, Szilassi's coordinates** [GS, Table 3], in exact rational arithmetic:
  - all seven faces are planar;
  - the exact checker says the polyhedron is embedded;
  - 8 of its 42 corners are reflex;
  - exactly one face, the one in the plane z = 4, is convex.

  Table 3 is not in general position: for example det(a₀, a₅, a₆) = 0 (one such zero for each pair of faces swapped
  by its half-turn symmetry). A perturbation of at most 5·10⁻⁴
  per plane coefficient makes it uniform while keeping it embedded and keeping all 42 reflex bits. On that
  configuration both encodings show no violated condition, and both encodings' reflex bits equal exact geometry at all
  42 corners.
- **C2 / D1:** with the 70 signs fixed to Szilassi's configuration, "at least 1 convex face" is SAT. "At least 2" and
  "none" are UNSAT, with verified proofs, in both encodings.
- **C3:** the signs of a random seven-plane configuration that the checker rejects give UNSAT. "At least 8 convex faces"
  is UNSAT. This is only a sanity check, since its cardinality clauses contain the empty clause.
- **C4 / D3:** both base instances are SAT. Demanding non-convexity of six named faces only is SAT, so the new clauses
  are not contradictory on their own.

**Oriented matroids (descriptive, `om_classes.py`).**
- Without any convexity clause, the ENC1 base instance has exactly 1,260 solutions (one per reorientation class),
  enumerated to exhaustion.
- Up to relabelling by the 42 automorphisms of the map and mirror image, they form 15 classes of 84 each. Szilassi's
  configuration is in one of them.
- Restricted to the seven face planes (dropping the plane at infinity), all 15 give the same oriented matroid. This
  agrees with the reported uniqueness of Szilassi's torus, and here it holds for all realisations, not only
  symmetric ones.
- The 15 classes differ in the position of the plane at infinity, which convexity depends on. The theorem shows
  each of them has exactly one convex face.
- Whether all 15 classes are realisable is not decided here.

**Independent review.** A referee re-derived the predicates by hand. It also tested about 10,000 embedded
configurations near Szilassi's (both handednesses) against its own exact embeddedness test, and found 0 violated
conditions and 0 reflex-bit mismatches. Its logs are not part of this bundle.

It also reported, post hoc and without a certificate here, that "at most one convex face" stays UNSAT after the
edge-piercing and one-reflex-corner-per-vertex clauses are removed. If so, the upper bound needs only face simplicity
and the plane arrangement, and a short hand proof may exist.

## 6. Reproduce

From this folder, with Python 3 and numpy, plus CaDiCaL and drat-trim on the path (or passed as `--cadical`,
`--drat-trim`):

```
python convex7.py checks          # C0, C1 (seconds)
python convex7.py cnf             # rebuild all 18 CNFs, compare SHA-256 with certificates.csv
python convex7.py prove           # solve + drat-trim every instance except the 4 slow ones (seconds)
python convex7.py prove --full    # all of them (about 10-15 more minutes, mostly drat-trim)
python om_classes.py              # the oriented-matroid count of §5 (needs python-sat)
```

The tools can also be set with the environment variables `CADICAL` and `DRAT_TRIM`. On Windows, both need the
binary-mode patch in `../PATCHES.md`. The checks use the exact checker `code/check.py` and the map from
`code/k12.py`.

## 7. Fairness note

This result agrees with Szilassi's coordinates [GS, Table 3]: one convex face, checked exactly in §5. It also agrees
with a sci.math post by Bill Taylor (8 May 1996) quoted in [EJ]: "the other one is convex and almost regular".

Bokowski and Schewe (2004 preprint, "On Szilassi's Torus", as reported in [SV]) studied the axially symmetric
realisations. The theorem here covers all realisations.

It disagrees with a sentence by Bokowski [BV, §1, paragraph above Figure 1]: "Szilassi's torus is shown on the right
with its seven non-convex facets." The theorem shows that no embedded seven-hexagon torus without overarching faces
has seven non-convex faces. The sentence most likely describes the figure loosely, and we do not claim it was meant
as a mathematical statement.

## References

- [GS] B. Grünbaum and L. Szilassi, Geometric realizations of special toroidal complexes, *Contributions to Discrete
  Mathematics* 4(1) (2009) 21–39, doi:10.55016/ojs/cdm.v4i1.61986. Definitions p. 21; Remarks 5.4–5.6 and Conjecture
  5.7 p. 32–33; Table 3 p. 37.
- [BV] J. Bokowski, On tori of Császár and Szilassi and more, VisMath electronic book,
  http://www.mi.sanu.ac.rs/vismath/visbook/bokowsky/, §1.
- [EJ] D. Eppstein, The Geometry Junkyard, "Szilassi polyhedron",
  https://www.ics.uci.edu/~eppstein/junkyard/szilassi.html (sci.math post by Bill Taylor, 8 May 1996, quoted there).
- [SV] L. Szilassi's page in the VisMath electronic book, https://www.mi.sanu.ac.rs/vismath/visbook/szilassi/index.html
  (on Bokowski and Schewe 2004).
- [AK] O. Arseneva, L. Kleist, B. Klemz, M. Löffler, A. Schulz, B. Vogtenhuber, A. Wolff, Adjacency graphs of
  polyhedral surfaces, *Discrete Comput. Geom.* 71 (2024) 1429–1455, arXiv:2103.09803, §3.1, Prop. 6.
- [BH] J. Bokowski and K. H., *Symmetry* 17(4) (2025) 622, doi:10.3390/sym17040622, §8.2.
- [MC] D. McCooey, Visual Polyhedra, Szilassi polyhedron versions 1 and 2,
  http://dmccooey.com/polyhedra/Szilassi1.html and http://dmccooey.com/polyhedra/Szilassi2.html ("Faces" row: 1 convex
  hexagon).
- drat-trim: https://github.com/marijnheule/drat-trim. CaDiCaL: https://github.com/arminbiere/cadical.
