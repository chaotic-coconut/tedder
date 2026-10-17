# tedder roadmap

Stage 00 — foundation repair and Stage 01 — kernels are complete. Their
detailed work is defined by the Stage 00
[field plan](docs/stages/00-foundation/field-plan.md) and [test plan](docs/stages/00-foundation/test-plan.md),
and the Stage 01 [kernel plan](docs/stages/01-kernels/kernel-plan.md).

Later stages are planned:

| Stage | Focus | Detailed document |
| --- | --- | --- |
| 02 — Reference neighbour search | Search contract and brute-force implementation | Not yet present |
| 03 — Degree 0 and 1 fitting | Fixed-bandwidth local fits | Not yet present |
| 04 — Degree 2 fitting | Quadratic fitting | Not yet present |
| 05 — nanoflann backend | kd-tree backend | Not yet present |
| 06 — Derived quantities | Divergence, vorticity, and strain | Not yet present |
| 07 — Optional I/O and policies | NumPy I/O and explicitly accepted policies | Not yet present |

Unscheduled candidates are not planned work:

| Candidate | When to reconsider | Notes |
| --- | --- | --- |
| Spherical scalar reconstruction | After Stage 06, and only for a concrete dataset with stated derivative-accuracy requirements | [ADR 0003, "Out of scope: spherical domains"](docs/adr/0003-geometry-local-frame-support.md#out-of-scope-spherical-domains) |
