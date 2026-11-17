# szilassi-12: figure captions and sources

All three figures are drawn by `figures/make_figures.py` (run: `python figures/make_figures.py .` from the repo root;
Python 3.11, numpy, scipy, matplotlib 3.11, a few seconds). Szilassi's coordinates are read as data from
`convex7/core.py` (`T3V`, `T3F`: Grünbaum & Szilassi 2009, Table 3) with the `ast` module, without importing or running
that module. Planarity, shared edges, inward (reflex) corners and convexity are recomputed in exact fractions.
Face numbers 1-7 follow the order of the faces in `T3F`.

## 1. `figures/szilassi-polyhedron.svg` (hero)

**Point.** Szilassi's polyhedron itself: a ring of seven flat hexagons in which every face touches every other, with
its one un-dented face picked out.

**Caption.** Szilassi's polyhedron, drawn from the published exact coordinates: seven flat hexagons, one colour each.
Every two of them share an edge (21 shared edges in all), and the solid is a ring with a hole you can see through
(front and back views). Six faces are dented, with one or two corners pointing inwards; face 7, the small face in the plane
z = 4 (seen in the view from below), is the one un-dented face. Part A of this repository proves that every
seven-hexagon ring of this kind with no two faces "overarching" (a technical condition, see Part A; rings with overarching faces are not covered) has exactly one un-dented face.

**Alt text.** Three views of a pointed, roughly tetrahedral solid made of seven coloured faces. The front and back
views each show a small see-through hole. The view from below shows a small lilac hexagon, labelled as the one
un-dented face. A legend lists the number of inward corners of each face: 1, 2, 1, 1, 2, 1 and 0.

**Sources.**
| element | file | field |
|---|---|---|
| 14 vertices, 7 faces | `convex7/core.py` | `T3V`, `T3F` (script asserts 14 and 7×6) |
| every pair of faces shares exactly one edge; 21 edges | computed from `T3F` | script asserts 21 edges, each pair once |
| inward corners per face 1, 2, 1, 1, 2, 1, 0 (8 of 42) | computed exactly (turn sign against face orientation, as in `convex7/core.py` `corners_exact`) | matches convex7/README §5 C1: "8 of its 42 corners are reflex" |
| face 7 in the plane z = 4 | computed (its normal is (0, 0, −115); all its vertices have z = 4) | matches convex7/README §5 C1 |
| views | computed | orthographic, one ray per pixel (900 px across), nearest face wins; az/el = (0°, 20°), (180°, 20°), (195°, −35°) |
| see-through hole | computed | background pixels not connected to the image border |

## 2. `figures/szilassi-faces-flat.svg`

**Point.** The theorem's statement made visible: six faces have inward corners, exactly one has none.

**Caption.** The seven faces of Szilassi's polyhedron laid flat, each in its own plane and all at the same scale.
Open circles mark corners that point inwards: 8 of the 42 corners. Faces 1 to 6 each have one or two; face 7 has none,
the polyhedron's one un-dented face. Part A proves that this count of un-dented faces, exactly one, is the same for
every seven-hexagon ring in which each face shares an edge with every other with no two faces "overarching" (a technical condition, see Part A; rings with overarching faces are not covered).

**Alt text.** Seven coloured hexagons in a row, drawn to the same scale. Six are irregular, with one or two white rings
marking inward corners; the seventh, much smaller, has none and is outlined in white.

**Sources.** As figure 1 (`T3V`, `T3F`; exact reflex test; counts 8 of 42).

## 3. `figures/szilassi-symmetry-cases.svg`

**Point.** Part B in one glance: every symmetric case of the 12-face question was proved impossible, each with a
re-checked proof file.

**Caption.** Part B: could a 12-face polyhedron in which every two faces touch have a symmetry? A symmetry would have to
be one of a short list of motions (a half-turn, a one-third turn, a quarter-turn combined with a mirror, or one of
the larger groups D2, D3 and T) acting on one of the 59 possible face patterns. Each dot is one such case, and its
position shows how long the solver took to prove it impossible (log scale): from 7 seconds to about 2 h 15 min. 25
cases (21 half-turns and 4 one-third turns) already settle the theorem; the larger groups were proved as extra checks.
All 36 proof files were re-checked by a separate checker program and came back VERIFIED. Whether a 12-face version
without any symmetry exists is still open.

**Alt text.** Dot chart with six rows: half-turn (21 dots, 19 of them between 40 minutes and about 2 h 15 min, two faster
ones at about 3 and 12 minutes), one-third turn (4), quarter-turn with mirror (6), D2 (2) and D3 (2), all between
about 25 seconds and 5 minutes, and T (1 dot at 7 seconds).

**Sources.**
| element | file | field |
|---|---|---|
| 36 cases, groups, times | `results/certificates.csv` | `group` (C2, C2#1, C2#2 counted as half-turn), `solve_s`, `verdict` (script asserts all VERIFIED) |
| 7 s and 8,114 s (about 2 h 15 min) | `results/certificates.csv` | min and max of `solve_s` |
| 25 cases = 21 half-turns + 4 three-fold | README "Part B in numbers"; TECHNICAL §2 | text |
| 59 face patterns | TECHNICAL §1 | text |
| the two fast half-turn cases (188 s, 735 s) | `results/certificates.csv` | maps 48 and 23, encoding `orbfix+orient` |

Dropped: a drawing of the face-adjacency graph (seven nodes, all joined). It is exact but adds nothing beyond the
sentence "every face touches every other", and the hero already states it.

Referee note (2026-11-17): the Part A sentences in captions 1 and 2 now carry the README's own scope ("no two of them
overarching"; "those *with* overarching faces are not covered", README At a glance, Part A). The hero caption in
PLACEMENT.md stops before the Part A sentence and needs no change.
