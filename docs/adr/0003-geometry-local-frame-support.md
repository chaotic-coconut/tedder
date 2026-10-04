# ADR 0003 — geometry, local frame and support

Status: proposed. Draft for the author's review; not authoritative until
accepted. On acceptance it supersedes the `distance_squared` requirement of
the `Domain` concept in [ADR 0001](0001-domain-periodic-contracts.md) F1.
Everything else in ADR 0001 and [ADR 0002](0002-point-type-restriction.md)
stands.

Needed before Stage 01 (kernels) and Stage 02 (neighbour search).

## Context

The implementation plan of 2026-10-04 and its
[review](../reviews/2026-10-04-implementation-plan-review.md) found four
geometric contracts that were undefined or could disagree:

- `Domain` defines both `offset` and `distance_squared`, and a custom model
  may make them disagree. `field-plan.md` F4 explicitly allowed custom models
  with another metric.
- Fitting needs local coordinates, but the plan never defined the support
  variable `q`, or which formula decides membership near `q = 1`.
- The search backend contract said only that it produces "possible hits".
- Anisotropic support was discussed, although nothing planned uses it.

## Decision

### 1. Displacement is the only geometric primitive

- `offset(p, q)` remains the one geometric operation a `Domain` defines. It
  returns `p − q`, reduced to the minimum image on periodic axes (ADR 0001 F3).
- `distance_squared` is removed from the `Domain` concept and from the
  built-in models.
- A free function `distance_squared(g, p, q)` returns the sum of the squared
  components of `g.offset(p, q)`, summed in axis order `0 … D−1`.
- `distance(g, p, q)` is the square root of that free function.
- For `Euclidean` and `Periodic` nothing changes numerically. Both member
  functions already compute the same sum in the same order.
- All other `Domain` requirements are unchanged.

| Part | Responsible for |
| --- | --- |
| `Domain` | Topology: how a displacement is obtained |
| Local frame (§3) | Scaling the displacement into local coordinates |
| Kernel | Influence as a function of the local distance |

tedder is Cartesian. A custom `Domain` may change how displacements are
formed, for example which periodic images exist. It cannot introduce another
metric. Manifold or spherical geometry needs its own ADR.

### 2. Bandwidth stays scalar and isotropic

`max_bandwidth()` and `admits_bandwidth(g, h)` keep their meaning and results
from ADR 0001 F2:

- `h` is admissible exactly when `0 < h < max_bandwidth()`;
- for `Periodic`, `max_bandwidth()` is half the shortest positive period;
- an invalid configuration admits no bandwidth.

Per-axis limits are not added to the public API. If a frame other than
`B = I/h` is introduced (§3), it is admissible only when, on every periodic
axis `d`:

```text
reach_d(B) = ‖row_d(B⁻¹)‖₂  <  L_d / 2
```

Here `reach_d` is how far the support `{r : ‖B r‖ < 1}` extends along axis
`d`. For `B = I/h` this reduces to `h < L_d / 2`.

No single scalar expresses the condition in general. An ellipse with
semi-axes 0.6 and 0.1, rotated 45°, reaches only ±0.43 along x. It is
therefore admissible for `L = 1`, even though one semi-axis exceeds `L / 2`.

The ADR that introduces such a frame must add the matching query. This ADR
only records the condition.

### 3. Local frame

For an evaluation point `x₀` and a sample `xᵢ`:

```text
rᵢ = g.offset(xᵢ, x₀)    physical displacement, wrapped
ξᵢ = B rᵢ                local coordinates
```

- **Order and inputs.** Wrapping happens before scaling. Sample coordinates
  are never modified.
- **Linearity.** `offset` performs the translation in `ξ = B(x − x₀)`, so `B`
  itself is linear.
- **Public frame.** Public behaviour uses only `B = I/h`.
- **Back to physical coordinates.** Derivatives return as below, with Jacobian
  rows indexed by component. For `B = I/h` this gives `J_x = J_ξ / h` and
  `H_x = H_ξ / h²`. ADR 0004 fixes how the Hessian transform keeps exact
  symmetry.

  ```text
  J_x     = J_ξ B
  H_x^(c) = Bᵀ H_ξ^(c) B      for each component c
  ```

- **Internal frame type.** The frame is a concrete type in namespace `detail`.
  Fitting code reaches the frame only through its operations: mapping an
  offset to `ξ`, and transforming the Jacobian and Hessians back. It never
  divides by `h` directly.
- **No public frame concept.**
  - One internal type handles isotropic scaling.
  - Tests may add a second internal type with a general matrix `B`.
  - Once both exist, a `detail`-only concept may constrain internal code.
  - A public frame concept needs its own ADR, together with anisotropic or
    space-time support.

### 4. Support

- **Support variable.** `q² = ξ · ξ`, computed from the canonical `ξ` of §3.
  For `B = I/h` the canonical form is `ξ_d = r_d / h`, one division per
  component.
- **Candidate test.** A sample is a candidate exactly when `q² < 1`. The
  comparison is written so that a NaN `q²` fails it.
- **Kernel input.**
  - The kernel receives `q` or `q²`, computed from the same `ξ`.
  - Where it needs `q`, `q = sqrt(q²)` is computed only for candidates.
  - For `float` and `double`, `q² < 1` holds exactly when `sqrt(q²) < 1`, so
    the two forms cannot disagree.
- **Kernel contract (Stage 01).** The kernel's weight is zero for `q ≥ 1`;
  Stage 01 defines the rest of the kernel contract.
- **Active samples.** A candidate whose weight is not finite and positive is
  not active (ADR 0004). Weight underflow near `q = 1` is one way this can
  happen.
- **No exact boundary.** Four candidate formulas were compared against exact
  arithmetic on 100,000 points near the boundary. Every one misclassified some
  points within about 3 ulps of `q = 1`.
- **Why this formula.** The canonical one is the kernel's own input, does not
  overflow or underflow for extreme `h`, and needs no square root for rejected
  points. "Authoritative" means every caller uses this one formula, not that
  it is exact.
- **Half-period ties.** In exact arithmetic, a sample exactly half a period
  away has `q ≥ L / (2h) > 1`. The tie convention of ADR 0001 F3 therefore
  never places an image inside the support.

### 5. Neighbour search

A search backend returns sample indices for one query. Its contract:

- **Complete.** It returns every sample with `q² < 1`, including samples
  reached across a periodic boundary.
- **Extras allowed.** It may return samples outside the support.
- **Duplicates allowed.** It may return an index more than once, for example
  when it queries several periodic images.
- **No decisions.** It does not decide candidacy and does not evaluate
  kernels.

The fitting layer then:

1. sorts the returned indices in ascending order and removes duplicates
   (`std::sort` and `std::unique` on a reused buffer);
2. applies the support test of §4;
3. counts the candidates.

Consequences:

- **Counts.** Candidate counts do not depend on the backend.
- **Row order.** Rows enter the design matrix in ascending sample index. For
  the same candidate set, every backend therefore gives the fitter identical
  input in identical order, and the Stage 05 equivalence tests compare fits
  with exact equality.
- **Search radius.** A backend may search with a slightly enlarged radius;
  the filter makes the final decision.
- **Periodic kd-tree (Stage 05).** The backend only has to be complete, for
  example by querying periodic images; the filter provides exactness.
- **Future frames.** For a frame other than `B = I/h`, a radius search needs
  the bounding radius `‖B⁻¹‖₂ = 1 / σ_min(B)`. This is recorded for the ADR
  that introduces such a frame.

## Consequences

- **`field.hpp`.**
  - Drop `distance_squared` from the `Domain` concept and remove the member
    functions.
  - Add the free `distance_squared`; `distance` calls it.
- **`tests/test_domain_contract.cpp`.**
  - Remove the fixtures `MissingDistanceSquared`, `DistanceSquaredWrongReturn`,
    `DistanceSquaredNotConst` and `DistanceSquaredNotNoexcept`. They would now
    satisfy `Domain`.
  - Remove the unused `distance_squared` member from the other fixtures, so
    each one still isolates a single defect.
- **`tests/test_domain.cpp`, `tests/test_geometry.cpp`.** Member calls become
  free-function calls. `test_geometry.cpp` already compares against an
  independently summed squared offset.
- **`README.md`.** Unaffected; it uses only `distance` and `max_bandwidth`.
- **On acceptance.** Update ADR 0001's status line, as was done for ADR 0002.
  Mark this ADR accepted in `docs/INDEX.md`, and link it from the Stage 01 and
  Stage 02 rows of `ROADMAP.md`.

## Rejected alternatives

- **Keeping `distance_squared` in the concept** with a documented requirement
  that it equal the squared offset. Agreement would become a promise the
  compiler cannot check.
- **Per-axis bandwidth limits in the public API now.** No support model uses
  them, and they would change ADR 0001 F2, the tests and the README for
  deferred functionality.
- **A public frame concept now.** There is only one model, so its
  requirements would be guessed.
- **`r · r < h²` as the support test.** It differs from the kernel input and
  fails for extreme `h`: once `h²` underflows to zero, even `r = 0` is
  rejected.
- **`std::set` for duplicate removal.** It allocates a node per hit and
  yields the same order as sort + unique.

## Open points

- Whether non-finite sample coordinates stay a precondition of fitting or
  become a checked failure. ADR 0004 covers non-finite values only.
- Whether the Stage 02 brute-force backend prefilters with an enlarged radius
  or returns every index.
