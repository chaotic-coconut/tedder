# ADR 0003 — geometry, local frame and support

Status: accepted on 2026-10-17, at the close of Stage 01. Supersedes the
`distance_squared` requirement of the `Domain` concept in
[ADR 0001](0001-domain-periodic-contracts.md) F1. Everything else in ADR 0001
and [ADR 0002](0002-point-type-restriction.md) stands.

Implementation: §1 is already on `main`, in commit e99f6a5. §2 describes the
current behaviour. The kernel contract of §4 was implemented in Stage 01
(kernels), following the
[Stage 01 kernel plan](../stages/01-kernels/kernel-plan.md). The rest of §3 to
§5 belongs to Stage 02 (neighbour search: `ξ`, the candidate test and the
search contract) and Stage 03 (fitting: transforming derivatives back).

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
- Library code calls it qualified, as `tedder::distance_squared`. An
  unqualified call would also find a function of the same name declared next
  to a user's domain type, and an exact non-template match would win over
  tedder's template, so the library would silently use the user's metric.
- For `Euclidean` and `Periodic` nothing changes numerically. Their former
  member functions computed the same sum in the same order.
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
  - The kernel receives `q = sqrt(q²)`, never `q²`. Wendland kernels contain
    odd powers of `q`, so they would need the square root anyway.
  - `q` is computed only for candidates.
  - `q²` is a sum of squares, so `q ≥ 0`. `sqrt` is correctly rounded and 1
    is representable, so `q² < 1` gives `q ≤ 1` in every rounding mode. The
    kernel's precondition therefore holds for every candidate.
  - In the default rounding mode, round to nearest, the bound is strict: for
    `float` and `double`, `q² < 1` holds exactly when `sqrt(q²) < 1`, so the
    two forms cannot disagree, and the fitter passes the kernel only `q` in
    `[0, 1)`. With upward rounding, `sqrt` of the largest value below 1
    rounds to 1. That candidate then gets `w(1) = 0` and is not active.
  - A NaN `q²` fails the candidate test and never reaches the kernel.
- **Kernel contract.** For a scalar type `T`, a kernel is an object `k` whose
  `const` call `k(q)` is `noexcept` and returns exactly `T`. The `Kernel`
  concept checks only this. A kernel that models the concept also promises,
  writing `w(q)` for `k(q)`:
  1. **Precondition.** `0 ≤ q ≤ 1`, which also excludes NaN. Any other `q`
     is a caller bug, not a reported error. The kernel may `assert` the
     precondition, as `Periodic::offset` asserts its own.
  2. **Boundary.** `w(1) = 0`. In round-to-nearest the fitter never passes
     `q = 1`. The closed interval lets tests check the boundary value, and
     keeps the precondition valid under directed rounding.
  3. **Range.** For `q` in `[0, 1]` the result is finite and `≥ 0`. It may
     round to 0 near `q = 1`, and the sample is then not active.
  4. **Positive inside.** In exact arithmetic, `w(q) > 0` for `0 ≤ q < 1`, so
     every candidate carries a weight unless that weight rounds to 0.
  5. **Scale.** Only relative weights matter. In exact arithmetic, multiplying
     every weight by the same positive constant changes neither the
     least-squares solution nor the weighted RMS of ADR 0004. A kernel
     therefore carries no normalising constant and does not depend on `h` or
     `D`.
  6. **Equality-preserving and free of side effects.** The result depends
     only on `q` and on state of the kernel that never changes, so the same
     `q` always gives the same result. A call writes nothing: not the kernel,
     not any object the kernel refers to, and no global state. Concurrent
     calls are therefore safe.

  The precondition matters because a kernel formula need not vanish, or even
  stay non-negative, outside `[0, 1]`. Wendland C2, `(1 − q)⁴ (4q + 1)`, gives
  `w(2) = 9`, nine times its value at the centre, and `w(−0.5) ≈ −5.06`.
- **Not in the kernel contract.**
  - The kernel's derivative. The fitter takes derivatives from the polynomial
    coefficients, never from `w′`.
  - Its smoothness order. It is a documented property of each kernel,
    verified on the kernel's formula in Stage 01.
  - The dimension `D`. Wendland's dimension parameter guarantees positive
    definiteness, which matters for RBF interpolation matrices.
    Least-squares weights only need to be non-negative.
  - `w(0) = 1`. The Stage 01 kernels follow this convention, but it is not
    required.
- **Active samples.** A candidate whose weight is not finite and positive is
  not active (ADR 0004). Weight underflow near `q = 1` is one way this can
  happen.
- **No exact boundary.** Four candidate formulas were compared against exact
  arithmetic on 100,000 points near the boundary. Every one misclassified some
  points within about 3 ulps of `q = 1`.
- **Why this formula.** The canonical one is the kernel's own input, and it
  needs no square root for rejected points. "Authoritative" means every
  caller uses this one formula, not that it is exact.
- **Range.** Near the boundary, every component of `ξ` has magnitude at most
  about 1, whatever the scale of `h`. So `ξ · ξ` cannot overflow there, and a
  component small enough to underflow is far too small to change the decision.
  Elsewhere the range is limited:
  - Far outside the support, `ξ · ξ` can overflow to infinity. The sample is
    still rejected, which is the right answer. Example: `float`, `h` the
    smallest normal value, `r = 1`.
  - Far inside the support, a component of `ξ` can underflow to zero, and the
    fit then treats that sample as lying at `x₀` along that axis. Example:
    `float`, `h` the largest finite value, `r` the smallest normal value,
    where `r / h` rounds to 0 instead of 3.5e-77.
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

## Out of scope: spherical domains

Spherical domains are an unscheduled candidate (see `ROADMAP.md`), not part
of this decision. A spherical domain would need its own ADR. These contracts
of this ADR, and the result types of ADR 0004, would have to be revised:

1. **Point dimension and fit dimension.** A point on the sphere is stored as
   3 coordinates (a unit vector), but the polynomial basis needs 2 independent
   local coordinates. `offset` could return the log map as a 3-component
   tangent vector, but those vectors all lie in one plane, so a polynomial
   basis of degree 1 or higher in all 3 coordinates gives a rank-deficient
   design matrix. The meaning of `D` throughout fitting and results would
   have to change, not only `point_type`.
2. **`offset` is no longer total.** Today `offset` is defined for every pair
   of points. The log map is not defined at the antipode, although the
   geodesic distance there is. A search backend, for example brute force, can
   return such a point before the support test rejects it. The §4 test would
   reject a NaN offset, but that is not a contract. A spherical design needs
   explicit out-of-chart behaviour, such as rejecting by geodesic distance
   before evaluating the log map. The support radius is a physical length,
   tied to the sphere's radius.
3. **What "physical coordinates" means.** For Cartesian domains, derivatives
   are returned in Cartesian components. On a sphere, a gradient is a tangent
   vector. It can be given as 2 components in a tangent basis, but no
   continuous tangent frame covers the whole sphere, so east/north breaks
   down at the poles. Or it can be given in ambient form, as a 3-component
   tangent gradient and a 3×3 Hessian acting on tangent directions. This is
   undecided, and today's result types do not promise to cover either form.
4. **Scalar channels and vector fields differ.** A map from the sphere to
   `R^C` can be fitted component by component. A tangent vector field has its
   samples in different tangent spaces and needs a transport rule before they
   can be combined. Its second covariant derivative is not symmetric in
   general: the two orders differ by the curvature acting on the vector. So
   the symmetric per-component Hessians of ADR 0004 cannot take on that
   meaning. Scalar fields would come first; tangent vectors and their first
   derivatives later; vector second derivatives would need their own
   definition.

What is already known: in geodesic normal coordinates centred at the
evaluation point, the connection coefficients vanish there. So the coordinate
gradient and Hessian of `f ∘ exp` at the centre are the intrinsic gradient and
covariant Hessian, and the trace of that Hessian is the Laplace–Beltrami
operator applied to `f` at that point. This was checked symbolically for
`f = xy + 0.3z` on the unit sphere. It is an identity at the centre, not a
statement about fitting accuracy at a finite bandwidth: `f ∘ exp` is
generally not a polynomial.

The sphere search in nanoflann 2.0 is a backend option only. It could serve
the completeness contract of §5, and it solves none of the issues above.

## Consequences

Done for §1 in commit e99f6a5:

- **`field.hpp`.** `distance_squared` is gone from the `Domain` concept and
  from `Euclidean` and `Periodic`. The free `distance_squared` sums the squared
  offset components, and `distance` calls it as `tedder::distance_squared`.
- **`tests/test_domain_contract.cpp`.**
  - The fixtures `MissingDistanceSquared`, `DistanceSquaredWrongReturn`,
    `DistanceSquaredNotConst` and `DistanceSquaredNotNoexcept` are removed,
    because they now satisfy `Domain`. The remaining fixtures no longer carry a
    `distance_squared` member.
  - Two tests guard the decision. "a member distance_squared is ignored by the
    free functions" uses a domain whose own member returns a wrong value.
    "distance ignores an ADL-visible distance_squared overload" puts a function
    of the same name next to a domain type. Making the call in `distance`
    unqualified makes the second test fail.
- **`tests/test_domain.cpp`, `tests/test_geometry.cpp`.** Member calls became
  free-function calls. `test_geometry.cpp` still compares against an
  independently summed squared offset, with exact equality.
- **`README.md`.** Unaffected; it uses only `distance` and `max_bandwidth`.

Done for the kernel contract of §4 in Stage 01:

- **`include/tedder/kernel.hpp`.** The `Kernel` concept, which checks the
  call with a `const T` argument, and `WendlandC2` and `WendlandC4`, evaluated
  in factored form with the precondition asserted.
- **`tests/test_kernel.cpp`.** Concept checks, and kernel values for `float`
  and `double` against exact reference weights from
  `tools/kernel_reference.py`.

Still to do:

- The rest of §3 to §5 in Stages 02 and 03. Their stage files will carry the
  tests.

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
- **`q²` as the kernel input.** It saves the square root only for kernels
  that are polynomials in `q²`, such as the biweight `(1 − q²)²`. Wendland
  kernels would have to take the root themselves.
- **A tagged kernel argument**, such as `Radius<T>`. It would turn a `q`
  versus `q²` mix-up into a compile error, but every call would need a
  wrapper, for a convention this ADR now fixes.
- **A total kernel** that returns 0 for `q > 1`. It repeats the support test,
  and a caller passing `q > 1` would go unnoticed instead of failing an
  assertion.

## Open points

- Whether to rename `max_bandwidth()`, for example to `bandwidth_limit()`. The
  value is a strict limit that is itself not admissible, and "max" suggests
  otherwise. It is optional; if done, do it before more code depends on the
  name.
- Whether non-finite sample coordinates stay a precondition of fitting or
  become a checked failure. ADR 0004 covers non-finite values only. Stage 02
  keeps them a precondition with a debug-only assertion
  ([Stage 02 plan](../stages/02-neighbour-search/search-plan.md), A5);
  Stage 03 decides.
- Decided in the Stage 02 plan (A1): the brute-force backend returns every
  index, without prefiltering. Candidate selection makes every decision.
