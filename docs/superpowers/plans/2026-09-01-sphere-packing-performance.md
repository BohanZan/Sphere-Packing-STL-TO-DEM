# SpherePacking performance execution record

This is the branch-local execution record for the performance plan. The
user-owned source plan exists outside this worktree and was not modified.

## Goal and invariants

Accelerate the STL-to-sphere-packing workflow while preserving the placement,
relaxation, refill, fixed-seed behaviour, output CSV schemas, and report
fields. Geometry shortcuts remain conservative: uncertain occupancy cells use
the exact parity test, and every candidate sphere still receives its finite
triangle clearance check.

## Task disposition

| Task | Status | Commit(s) |
| --- | --- | --- |
| 1 — baseline fixtures and benchmark | completed | `3ee467c` |
| 2 — collision-before-ray reorder | reverted by user instruction | `b3358b1` was reset; not part of this branch |
| 3 — cached inward normals/refill guard | completed | `6061955` |
| 4 — batched exact vertical-ray parity | completed | `04d5596` |
| 5 — conservative occupancy triage | completed | `3aa0510` |
| 6 — batched finite-triangle distance checks | completed | `f169bc8` |
| 7 — sphere-neighbour allocation reduction | completed | `2d82973` |
| 8 — preallocated state and incremental reindex | completed | `895d3eb` |
| 9 — final evidence and documentation | completed | pending documentation commit |

## Final evidence (2026-09-02)

The unprofiled three-run iron-particle benchmark used the saved fixed
configuration and seed 53. The retained pre-optimization baseline was
`126.735856, 128.264100, 146.996799` seconds (median `128.264100`), with 20
accepted spheres in all runs. Final measurements were `72.507578, 73.310399,
75.290296` seconds (median `73.310399`), again with 20 accepted spheres and
sphere-volume fraction `-0.237247632907` in every run.

The final median is 42.84% lower than the retained baseline. The baseline
artifact did not preserve CPU/MATLAB metadata, so this is recorded as the
saved fixed-case comparison rather than a cross-environment guarantee. Final
environment: 13th Gen Intel(R) Core(TM) i7-13850HX; MATLAB R2025a Update 1,
version 25.1.0.2973910.

The post-documentation full regression suite passed 26/26 with no failures or
incompletes in 11.4865 seconds. A capped four-sphere profile recorded 19,918 calls;
`verticalTriangleHit` and `spPointTriangleDistance` recorded zero calls.
The remaining targeted totals were 391 for `spHashNeighbours`, 391 for
`spSphereHitsTriangles`, and 258 for `spExactPointInside`.

The one-run profiler comparison of the iron-particle case showed that enabling
occupancy triage preserved accepted count and fraction while reducing exact
point-inside calls from 60,745 to 29,421 (51.57%). Its profiler-on wall time is
not used for the benchmark claim.

## Residual work

The reduced profile still shows CSV writing and relaxation work. The retained
triangle/hash/exact-call totals above are diagnostic evidence only; another
optimization should begin with a representative profile and a fresh
equivalence test, not by changing packing logic.
