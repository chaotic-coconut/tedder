# Review: proposed implementation plan of 2026-10-04

Status: review notes, not authoritative. Nothing here is accepted until it is
recorded in an ADR. Written in reviewer mode; no production code was changed.

Outcome: the decisions are carried forward into the draft
[ADR 0003](../adr/0003-geometry-local-frame-support.md) and
[ADR 0004](../adr/0004-fit-results-diagnostics-errors.md). Where those ADRs
differ from these notes, the ADRs take precedence. Changes made after the
review:

- D5: per-axis bandwidth limits are deferred. The condition is only
  documented.
- D14: non-finite sample values are rejected by default instead of dropped.
- D13: the return type is `std::expected<FitSuccess, FitFailure>` (shape B).
- D3: the SVD-of-R value is an estimate of κ₂, not an exact value.
- D15: the weighted RMS is a fit-quality measure, not a noise estimate.
- G8: rank invariance holds only for a fixed active set with positive weights.
- G6: the degeneracy tests cover any conic, not only circles.

Reviewed document: *tedder implementation plan after architecture discussion*
(2026-10-04, status "proposed"). It is not in the repository; `§` numbers below
refer to its sections.

Read for this review: `AGENTS.md`, `docs/INDEX.md`, `README.md`, `ROADMAP.md`,
`docs/AI_DEVELOPMENT.md`, ADR 0001, ADR 0002, both Stage 00 plans,
`include/tedder/field.hpp`, `tests/test_types.cpp`. The checks behind each
numerical claim are listed in [section 8](#8-checks-run).

## 1. Summary

The plan's direction is sound and compatible with the existing contracts.

- The Jacobian layout `jacobian[c, d] = ∂F_c/∂x_d` and the sign of `offset`
  are unchanged.
- ADR 0001's rescaling rule `J_x[:, d] = J_u[:, d] / s_d` is the diagonal case
  of the plan's `J_x = J_ξ B`.
- With the strict limit `h < L/2` and the strict support law `q < 1`, a
  half-period tie lies outside the support, so ADR 0001's tie-sign convention
  cannot affect a fit.

The formulas are correct: the gradient and Hessian transforms, the
streamfunction Jacobian, `ω = −Δψ`, `W = −4 det H_ψ`, and
`κ₂(AᵀA) = κ₂(A)²`.

The problems are elsewhere:

- some claims are stronger than floating-point arithmetic supports;
- several contracts are left undefined;
- the plan does not yet map onto the repository's documents and stages.

Most of these were settled in discussion. [Section 3](#3-decisions) records the
decisions and [section 5](#5-open-decisions) lists what remains open.

## 2. Findings on the plan

### 2.1 Fit with the repository

| ID | Finding | Resolution |
| --- | --- | --- |
| R1 | §25 asks to update "the authoritative project design document", but `PROJECT.md` does not exist. It is referenced by `docs/AI_DEVELOPMENT.md:4`, `test-plan.md:177` and `field-plan.md:201`. | Split the plan into ADRs, `PROJECT.md` and stage files; see [section 7](#7-documentation-changes). |
| R2 | The phases skip Stage 01 (kernels) and Stage 02 (search). The support law, the definition of `q`, and the candidate counts are kernel and search contracts. Phase B is Stage 03 work, and Phase D is not on the roadmap. | Assign each decision to the first stage that needs it; see [section 6](#6-placement-by-stage). |
| R3 | The existing `LocalFit` (`field.hpp:231`) has `neighbors` and `bandwidth`, and `tests/test_types.cpp:180-181` tests both. §3 forbids counts in `LocalFit`, and a scalar bandwidth stops being meaningful once `B` is general. | `LocalFit` holds only field information. Counts move to `FitDiagnostics` (D12). The bandwidth is an input the caller already has. |
| R4 | The linear-algebra dependency is not mentioned, although `AI_DEVELOPMENT.md:394` asks Stage 03 to record it. | Eigen (D2). |

### 2.2 Mathematics and numerics

| ID | Finding | Resolution |
| --- | --- | --- |
| N1 | §5.3 promises exact Hessian symmetry. `Bᵀ H_ξ B` sums the `(i, j)` and `(j, i)` entries in different orders, so they can differ by a few ulps. In 9,330 of 10,000 random 3×3 cases the result was not bitwise symmetric. | D1 |
| N2 | §9.2 says the ratio of R's diagonal entries is "not exact κ₂". More precisely it is a **lower bound** and can be very loose: for the 30×30 Kahan matrix the ratio is 3.42 while κ₂ is 9.88e3. | D3 |
| N3 | "Wrap first, then transform" (§5.1) is only safe while the support fits inside half a period on each periodic axis. `h < L/2` ensures this for isotropic scaling; a general `B` needs a per-axis condition. | D4, D5 |
| N4 | `q` is never defined. `Domain` allows `distance_squared ≠ ‖offset‖²` (`field-plan.md:188`). The candidate formulas also round differently near `q = 1`. | D6, D8 |
| N5 | §13 says the backend produces "possible hits". The contract also needs completeness and a rule for duplicate indices. | D9 |
| N6 | `B` is both the coordinate transform (§5) and the coefficient matrix in `AB ≈ Y` (§6). §4 writes `p(r)`, but the basis is built in `ξ`; the coefficients are ξ-derivatives. | Rename one `B`, and state that coefficients are ξ-derivatives. |
| N7 | §7 claims derivative observations would leave coefficient extraction and results unchanged. In the streamfunction case, each sample adds two rows for one unknown, the value is not identifiable, and derivative rows need `B` and carry other units. | D10 |
| N8 | "Order" (plan) and "degree" (roadmap) are used for different things, and it is unclear whether `LocalFit`'s parameter is the polynomial degree or the highest returned derivative. | D11 |
| N9 | The Taylor ½ factor has no effect for degree ≤ 1, so it is not urgent for Stage 03. Stage 03 does need a basis ordering, which the plan asks for without proposing one. | Suggested ordering: total degree first, then lexicographic `i ≤ j`. Basis size is `C(D+n, n)`: 1, `D+1`, `(D+1)(D+2)/2`. |

### 2.3 Gaps

| ID | Gap | Resolution |
| --- | --- | --- |
| G1 | `FitFailure::numerical_rank` is a plain `size_t`, but no factorization exists for validation failures or too few samples. | D12 |
| G2 | Non-finite sample values have no rule, and all components share one design matrix. | D14 |
| G3 | The residual norm depends on the kernel's normalising constant. | D15 |
| G4 | `float` is tested today, but the ε sweep down to 1e-12 and fixed thresholds only suit `double`. | D16 |
| G5 | A successful fit has no channel for diagnostics. | D12, D13 |
| G6 | The Stage 04 tests miss the realistic quadratic degeneracy: points on a conic. Twelve points on a circle give rank 5 of 6, because `x² + y²` is constant there. Ring stencils and polar grids produce this. | Add circle (2D) and sphere (3D) layouts to the Stage 04 tests. |
| G7 | §16: the Q criterion needs its norm stated (Frobenius). The Laplacian needs Hessians, so Stage 06 depends on Stage 04, while the roadmap says "from value/Jacobian data". | State the norm. Use `Q = −W/4`, valid for 2D divergence-free flow, as an independent test oracle. Update the roadmap. |
| G8 | Not recorded: an invertible `B` never changes the exact rank, because polynomials of total degree ≤ n are closed under affine maps. Full rank is therefore a property of the sample layout. Numerical rank does depend on `B`, which is the real reason to scale `ξ` to order 1. | Record in the fitting ADR. |

## 3. Decisions

"Decided" means the author decided it in discussion. "Recommended" means it is
a suggestion that still needs the author's choice. Both still need an ADR.

| ID | Topic | Status |
| --- | --- | --- |
| D1 | Hessian symmetry: compute `i ≤ j`, mirror | Recommended |
| D2 | Eigen, kept internal and replaceable | Eigen decided; the rest recommended |
| D3 | Exact κ₂ from an SVD of R | Recommended |
| D4 | The local transform applies only to local offsets | Decided |
| D5 | Per-axis bandwidth limits | Decided; shape recommended |
| D6 | Geometry derived from `offset`; `distance_squared` leaves the concept | Agreement decided; mechanism recommended |
| D7 | No public concept for transforms yet | Recommended |
| D8 | Support test `q² = ξ·ξ < 1` | Recommended |
| D9 | Complete search; duplicates removed by sort + unique | Removal decided; mechanism recommended |
| D10 | Keep the observation-row seam, drop the overclaim | Decided |
| D11 | `LocalFit`'s parameter is the highest returned derivative | Decided |
| D12 | Diagnostics in a separate, stage-wise struct | Decided |
| D13 | Use `std::expected` | Decided; toolchain open |
| D14 | Non-finite values: drop the sample and report it | Decided |
| D15 | Weighted RMS plus degrees of freedom | Recommended |
| D16 | `float` covered throughout | Decided |

### D1. Hessian symmetry

The `(i, j)` and `(j, i)` results are two roundings of the same exact number;
neither is more accurate. Use one rule: compute the entries with `i ≤ j` and
mirror them.

- It is cheaper: `D(D+1)/2` entries instead of `D²`.
- It is deterministic, and the stored matrix is exactly symmetric, so tests can
  compare with `==`.
- Averaging also gives symmetry, but costs more and gains no accuracy.

Isotropic scaling multiplies every entry by the same `1/h²`, so symmetry
survives without the rule. It only matters once `B` is general. In Eigen,
`selfadjointView<Upper>()` reads the upper triangle as the whole matrix.

### D2. Eigen without locking out LAPACK, MKL or Blaze

- **No Eigen in public signatures.** Public types stay tedder's own (`Point`,
  `Matrix`, `LocalFit`, diagnostics). Eigen lives in `detail`, behind one
  narrow internal step: factor, solve the right-hand sides, report rank.
- **`field.hpp` stays free of Eigen.** Fitting goes in its own header, so code
  that only uses geometry does not need Eigen.
- **MKL and LAPACK routing belong to the consumer.** `EIGEN_USE_MKL_ALL` and
  `EIGEN_USE_LAPACKE` send some decompositions to those libraries (for example
  `ColPivHouseholderQR` to `?geqp3`). They change Eigen for the whole
  translation unit, so a header-only library must never define them.
- **Column-major design matrix.** Householder QR works on columns, and LAPACK
  expects column-major data. A later LAPACK backend could pass the buffer to
  `?gelsy` (pivoted-QR least squares) or `?gelsd` (SVD least squares).
- **Rank tolerance in tedder's own terms.** For example `|r_kk| ≤ τ·|r_11|`,
  where τ is a documented multiple of `ε(T)`. If it were "Eigen's default",
  the contract would change with the backend.
- **No backend abstraction until a second backend exists.** This follows
  `AI_DEVELOPMENT.md`.
- **Benchmark before expecting MKL to help.** These matrices have at most 10
  columns and tens of rows, and at that size Eigen's inline code often beats a
  library call.

### D3. Condition numbers

For a full-rank `A` with `m` columns and `AP = QR`:

```text
|r₁₁ / r_mm|  ≤  κ₂(A) = σ_max/σ_min  ≤  ‖R‖_F · ‖R⁻¹‖_F  ≤  m · κ₂(A)
```

- The SVD ratio is the true value; the other two are bounds.
- For a warning, an underestimate is the dangerous direction: it reports
  "fine" when the problem is not. An upper bound can only raise false alarms.
- `A` and `R` have the same singular values, because `Q` has orthonormal
  columns and `P` is a permutation. `R` is at most 10×10, so an SVD of `R`
  after the QR gives the exact κ₂ for a few thousand flops. Eigen's
  `JacobiSVD` uses the same QR preconditioning internally for tall matrices.
- Exact κ₂ can therefore be an optional diagnostic. A full SVD of `A` is only
  needed as a test oracle.
- Each condition value carries its provenance, because a pivoted-QR bound, an
  SVD value and a normal-equations estimate are not equally reliable.

### D4. The transform applies only to local offsets

Sample coordinates never change. `Domain::offset` wraps first; `ξ = B·r` is
applied only to that local displacement. The wrap is always safe. The risk is
in which periodic images fall inside the support.

Example: x is periodic with `L = 1` and `B = diag(1/0.7, 1/0.1)`, so the
support reaches ±0.7 along x. A sample at `r_x = 0.4` has a second image at
`r_x = −0.6`, which is also inside. The minimum image keeps only the first, so
the fit no longer sees the neighbourhood the support describes.

With isotropic scaling this cannot happen, because `h < L/2`. A non-isotropic
`B` needs the per-axis condition in D5.

Note that the translation part of `ξ = B(x − x₀)` is already done by
`Domain::offset`, so what remains is a linear map, not an affine one.

### D5. Per-axis bandwidth limits

Replace the scalar `max_bandwidth()` in the `Domain` concept with a per-axis
query. The value for axis `d` is a strict upper limit on how far the support
may reach along that axis.

| Domain | Limit on axis `d` |
| --- | --- |
| `Euclidean` | `+∞` |
| `Periodic`, valid, `L_d > 0` | `L_d / 2` |
| `Periodic`, valid, `L_d ≤ 0` (open axis) | `+∞` |
| `Periodic`, invalid | `0` on every axis |

Recommended shape: one `noexcept` member returning all `D` limits at once, as
`point_type` (that is, `std::array<T, D>`). An indexed form
`bandwidth_limit(std::size_t d)` would need a `d < D` precondition, and every
consumer needs all axes anyway.

The current scalar becomes a free function, the minimum over the axes.
`admits_bandwidth(g, h)` keeps its current results. The values are identical
to today's `max_bandwidth()`: the per-axis values are the same `L_d / 2`
computations, and taking a minimum is exact.

Returning zero on every axis for an invalid domain keeps ADR 0001 F2's rule
that an invalid configuration admits nothing, for any future support shape as
well as the isotropic one.

Future anisotropic admissibility, not to be added now:

```text
reach_d(B) = ‖row_d(B⁻¹)‖₂  <  limit_d   for every axis d
```

No single scalar can express this. Example: an ellipse with semi-axes 0.6 and
0.1, rotated 45°, reaches only ±0.43 along x. It is admissible for `L = 1`
even though one semi-axis exceeds `L/2`.

Impact:

- public API: `max_bandwidth` appears on 4 lines of `field.hpp`, 1 line of
  `README.md`, and 56 lines across the tests;
- documentation: supersedes the `max_bandwidth` wording of ADR 0001 F2.

### D6. Geometry derived from `offset`

To make agreement impossible to break rather than merely documented:

- Remove `distance_squared` from the `Domain` concept.
- Provide a free `distance_squared(g, p, q)` that sums the squared components
  of `g.offset(p, q)` in axis order.
- `distance` uses that free function.

Both built-in domains already compute exactly this sum in this order, so
results are bitwise unchanged.

This also separates responsibilities, in line with the plan's §19:

- `Domain` describes topology: how to obtain a displacement.
- The local frame (`B`) describes scaling and the support shape.

Custom domains can no longer define a non-Euclidean metric. That is intended,
since tedder stays Cartesian.

Impact: this supersedes part of the concept accepted in ADR 0001 F1.
`distance_squared` appears on 6 lines of `field.hpp` and 35 lines of the tests.

### D7. No public concept for transforms yet

Recommendation: do not add one now; keep an internal seam instead.

The operations the fitter needs from a local frame are what a concept would
eventually list:

1. map a local offset to `ξ`;
2. convert the Jacobian back: `J_x = J_ξ B`;
3. convert the Hessians back: `H_x = Bᵀ H_ξ B`, using the D1 rule;
4. per-axis reach `‖row_d(B⁻¹)‖₂`, for admissibility (D5);
5. a bounding radius `‖B⁻¹‖₂ = 1/σ_min(B)`, for an isotropic search backend
   that must return a superset (D9).

Why not a concept now:

- There is one model (isotropic scaling), so the requirements would be
  guessed.
- A public concept is an extension point. `AGENTS.md` says not to turn
  undecided policies into public APIs, and `AI_DEVELOPMENT.md` says to
  generalise only when a second model creates pressure.

What to do now:

- Write one concrete internal type in `detail` (isotropic scale) that provides
  the operations above.
- The fitter calls only those operations and never divides by `h` directly.
- The plan wants a non-diagonal transform in internal tests (§23, Phase E).
  A second internal type (a general matrix) used only by tests provides it.
- Once those two exist, a `detail` concept constraining the internal fitter
  is reasonable, because it stays private.

Make it public only when anisotropic or space-time support is accepted by its
own ADR. Naming hint: the object fixes the metric, the scaling and the
admissible reach together, so "local frame" or "support" may describe it
better than "affine transform".

### D8. Support test

I compared four formulas against exact rational arithmetic, on 100,000 points
within `2e-15` of the boundary (double):

| Formula | Misclassified | Farthest misclassified, in units of 2⁻⁵³ in `q²` |
| --- | ---: | ---: |
| `Σ (r_d / h)² < 1` | 1051 | 2.0 |
| `r·r < h·h` | 1120 | 2.3 |
| `Σ (r_d · (1/h))² < 1` | 1440 | 3.2 |
| `sqrt(r·r) / h < 1` | 2024 | 3.5 |

Every formula is wrong only within about 3 ulps of the boundary, so accuracy
does not decide. Use `q² = ξ·ξ < 1`, with `ξ = B·r` (for isotropic scaling,
`r_d / h`):

- **It is the kernel's own input**, so a point cannot be a candidate by one
  formula and outside the support by another.
- **It does not overflow or underflow at extreme `h`.** `ξ` is of order 1 near
  the support. By contrast, `r·r < h²` rejects even `r = 0` once `h²`
  underflows to zero.
- **Rejected points skip the square root.** Taking the root after the test is
  safe: `q² < 1` holds exactly when `sqrt(q²) < 1`, for both `float` and
  `double`, because the square root of the largest value below 1 stays below 1.
- **Write it as `q2 < 1`, not `!(q2 >= 1)`**, so that a NaN offset fails the
  test.

Use one canonical function for `ξ`. The search backend may use a slightly
larger radius; this filter makes the final decision.

### D9. Search completeness and duplicates

Stage 02 contract:

- the backend must never miss a point in the support, including across a
  periodic seam;
- it may return extra points;
- it may return the same index more than once, for example when it queries
  several periodic images.

The filter removes duplicates before anything is counted. Collect hit indices
in a reused `std::vector<std::size_t>`, then `std::sort` + `std::unique` +
erase.

- **Why not `std::set`:** it allocates one node per hit and ends up in the
  same sorted order.
- **Backend-independent row order:** rows enter the design matrix in ascending
  sample index whatever the backend. The same candidate set then gives a
  bitwise identical fit, so the Stage 05 brute-force-versus-kd-tree test can
  use exact equality.
- **Alternative for huge clouds with many queries:** a per-sample stamp array,
  bumped once per query, avoids the sort at the cost of `N` integers per
  thread. Start with sort + unique.

This contract also settles the open kd-tree/periodic question in
`AI_DEVELOPMENT.md` for Stage 05: the backend only has to be complete.

### D10. The observation-row seam

Keep it and drop the claim that it leaves extraction and results unchanged.
Four rules keep future derivative observations possible without building them:

- row construction is separate from factorization;
- the solver does not assume one row per sample;
- the solver accepts rank-deficient systems and reports the rank, and the
  fitting policy decides what is acceptable;
- the internal steps produce coefficients plus a solver report, and
  `LocalFit` is assembled from them in a separate step. A future gauge-free
  fit could assemble a different result type from the same coefficients.

### D11. Polynomial degree and derivative order

There are two separate numbers:

- **degree `n`**, the polynomial that is fitted;
- **derivative order `K`**, what is returned: 0 is the value, 1 adds the
  Jacobian, 2 adds the Hessians.

`K ≤ n` is required. The natural default is `n = K`, but `n = K + 1` often
gives a better K-th derivative. Example: `f(x) = x²/2`, whose true slope at
`x₀ = 0` is 0.

| Sample points | Linear-fit slope | Quadratic-fit slope |
| --- | ---: | ---: |
| One-sided `[0, 0.1]` | +0.05 | ≈ 0 |
| Symmetric `[−0.1, 0.1]` | ≈ 0 | ≈ 0 |

On a one-sided neighbourhood, the curvature the line cannot represent leaks
into the slope, with a bias of about `c·h/2`. Every point near the edge of the
data has such a neighbourhood. The quadratic fit takes the curvature into its
own coefficient, at the cost of more unknowns: it needs more points and is
noisier. Local-polynomial theory (Fan & Gijbels, 1996) prefers an odd
`n − K` for this reason.

Decision: `LocalFit<T, D, C, K>` is indexed by the derivative order `K`. The
fitting function takes the degree separately, constrained by `degree ≥ K`.
Start by allowing only `degree = K`; allowing `K + 1` later changes no result
type. Name the parameter after what it means (for example `DerivativeOrder`),
and reserve "degree" for the polynomial.

### D12. Diagnostics in a separate struct

What exists at the moment of failure determines which fields can be present:

| Stage | Possible failures | Available information |
| --- | --- | --- |
| Validation | invalid domain, rejected bandwidth, non-finite `x₀` | none |
| Search and filter | too few active samples | counts |
| Factorization | rank deficient | counts, solver report |
| Solve | non-finite coefficients | everything |

So diagnostics fill up stage by stage, with one optional member per stage.
Conceptual sketch; names and types are undecided:

```cpp
struct FitCounts {
  std::size_t candidates;        // q² < 1, after duplicate removal
  std::size_t nonfinite_dropped; // candidates with a non-finite value (D14)
  std::size_t zero_weight;       // remaining candidates whose weight is 0
  std::size_t active;            // rows of the design matrix
};

template <Real T> struct SolverReport {
  std::size_t numerical_rank;
  std::optional<T> condition_estimate;
  ConditionEvidence evidence;    // none, cpqr_lower_bound, svd_exact, ...
};

template <Real T, std::size_t C> struct FitDiagnostics {
  std::optional<FitCounts> counts;              // once the search ran
  std::optional<SolverReport<T>> solver;        // once a factorization ran
  std::optional<std::array<T, C>> weighted_rms; // once a solve ran (D15)
};
```

- On success, every member is present; document this as an invariant.
- Do not store what can be derived. `basis_size` is a compile-time function of
  `D` and the degree, so expose it as `constexpr`. Degrees of freedom equal
  `active − numerical_rank`.

### D13. `std::expected`

Possible error codes: invalid domain, rejected bandwidth, non-finite
evaluation point, insufficient samples, rank deficient, non-finite result.

Two possible return shapes:

```cpp
// (A) Diagnostics beside the result: present whatever the outcome.
struct FitOutcome {
  std::expected<LocalFit<T, D, C, K>, FitError> result;
  FitDiagnostics<T, C> diagnostics;
};

// (B) Diagnostics inside both arms.
std::expected<FitSuccess, FitFailure>;
// FitSuccess{ fit, diagnostics }, FitFailure{ code, diagnostics }
```

- **(A)** keeps the error type small and makes the diagnostics reachable
  without branching. The cost is that the function returns a struct, so
  `and_then` chains start from `.result`.
- **(B)** returns `std::expected` directly, but stores the diagnostics in both
  arms.

Either way, D14's "dropped samples but still a result" case returns its
diagnostics. Recommendation: (A).

Toolchain consequence: see [section 4](#4-stdexpected-on-the-supported-toolchains).

### D14. Non-finite values

All `C` components share one design matrix (geometry); each component is one
right-hand-side column. Take sample 7 with value `(1.0, NaN)` and `C = 2`:

| Option | Consequence |
| --- | --- |
| Drop it for component 1 only | The components get different design matrices and need separate factorizations, so the shared design is lost. |
| **Drop it for all components** | One design matrix and one factorization; component 0 loses one good observation. |
| Do nothing | The NaN reaches only column 1 (each column is solved independently). Detectable, but silent. |
| Fail | No result. |

Decision: drop the sample for all components and report it in
`FitCounts::nonfinite_dropped`. The counts then read:

```text
active = candidates − nonfinite_dropped − zero_weight
```

Dropping can push a fit into "too few samples" or rank deficiency, which then
fails normally with the drop visible in the diagnostics. Whether any drop
should count as an error remains the caller's policy.

Non-finite coordinates are still a precondition of `Periodic::offset`. The
`q2 < 1` test (D8) would exclude a NaN offset automatically, but relying on
that means relaxing the precondition and counting those drops too.

### D15. Residuals

Weighted least squares minimises `Σ wᵢ eᵢ²`. Multiplying every weight by a
constant `c` does not change the solution, but the residual norm
`sqrt(Σ wᵢ eᵢ²)` grows by `√c`. Kernels often carry a normalising constant
that depends on `h` and `D`. The same data and the same fit would then report
different norms just because of the kernel convention or the bandwidth.

Example with `e = (0.1, −0.2, 0.1)` and `w = (1, 0.5, 0.25)`:

| Weights | `‖r‖` | Weighted RMS |
| --- | ---: | ---: |
| `w` | 0.1803 | 0.1363 |
| `100·w` | 1.8028 | 0.1363 |

Report, per component, the weighted RMS `sqrt(Σ w e² / Σ w)`. It is a weighted
average error, in the units of the values, and comparable with measurement
noise. Also report the degrees of freedom: when they are zero, the residual is
zero by construction and says nothing.

### D16. `float`

- Every tolerance and rank threshold is a multiple of
  `std::numeric_limits<T>::epsilon()`, never a literal such as `1e-12`.
- The degeneracy sweep runs down to about `ε(T)`: 1e-1 to 1e-7 for `float`.
- Accuracy runs out early. With κ₂ around 1e4, `float` Hessians may keep only
  about 3 correct digits, so reproduction tolerances should scale with
  `κ₂ · ε(T)`.
- `1/h²` overflows `float` when `h` is below about `5e-20`. Document the range
  rather than guarding it.
- Keep an internal `compute_type` alias, equal to `T` for now. If `float`
  proves too inaccurate, computing in `double` for `float` input becomes a
  one-line change.

## 4. `std::expected` on the supported toolchains

libstdc++ only provides `<expected>` when `__cpp_concepts >= 202002L`, and
Clang 18 defines `201907L`. Tested in this container; Ubuntu 24.04 packages:

| Toolchain | `__cpp_concepts` | `std::expected` | tedder suite |
| --- | --- | --- | --- |
| GCC 13.3, libstdc++ 13 | n/a | available | not run here; CI covers it |
| Clang 18.1.3, libstdc++ 13 (current CI) | `201907L` | **not available** | n/a |
| Clang 18.1.3, libc++ 18 | `201907L` | available | 50/50 passed |
| Clang 19.1.1, libstdc++ 13 | `202002` | available | 50/50 passed |

Options:

1. **Raise the minimum Clang to 19.** CI installs `clang-19`; both compilers
   keep libstdc++. This is the simplest option, and the README requirements
   change.
2. **Keep Clang 18, but with libc++.** CI installs `libc++-18-dev` and builds
   with `-stdlib=libc++`. This adds a second standard library to the test
   matrix. Consumers using Clang 18 with libstdc++ still cannot use the
   fitting header.

Recommendation: option 1, optionally with a libc++ job as extra coverage. In
the fitting header, check `__cpp_lib_expected` (from `<version>`) and emit a
readable `#error`, rather than leaving users with "no template named
'expected'". `field.hpp` does not need `std::expected` and can keep supporting
Clang 18.

## 5. Open decisions

- Return shape (A) or (B) (D13), and the Clang option (section 4).
- Names: the per-axis limits member, the free scalar function, the derivative
  order parameter, and whether the fourth `LocalFit` parameter has a default.
- Confirm the basis ordering (N9).
- The value of the rank tolerance τ, from the degeneracy experiments in §22,
  separately for `float` and `double`.
- Whether `compute_type` stays `T` for `float` (D16).
- How Eigen is found: `find_package(Eigen3 CONFIG)`, `FetchContent`, or both;
  and the minimum Eigen version.
- Whether non-finite coordinates move from precondition to counted drop (D14).

## 6. Placement by stage

| Before | Decisions |
| --- | --- |
| Stage 01, kernels | D4, D5, D6, D8: the kernel receives `q` from `ξ = B·offset`, `w(q) = 0` for `q ≥ 1`, per-axis limits, geometry derived from `offset` |
| Stage 02, search | D9: completeness and duplicate removal; candidate counts |
| Stage 03, degree 0 and 1 | R3, D2, D3, D7, D10, D11, D12, D13 with the toolchain choice, D14, D15, D16; basis ordering (N9) |
| Stage 04, degree 2 | D1, Taylor ½ factor, conic degeneracy tests (G6) |
| Stage 06, derived quantities | G7: Frobenius norm in Q, `Q = −W/4` oracle, dependency on Stage 04 |

## 7. Documentation changes

- **ADR 0003, geometry and support** (D4, D5, D6, D8, D9). Supersedes the
  `distance_squared` requirement of ADR 0001 F1 and the `max_bandwidth`
  wording of F2. Accept before Stage 01.
- **ADR 0004, fit results, diagnostics and errors** (R3, D10–D16, N9,
  toolchain). Accept before Stage 03.
- **ADR 0005, linear algebra** (D2, D3, the rank tolerance definition). This
  could be merged into 0004. Accept before Stage 03.
- **`PROJECT.md`**, new: the explanatory material of §16–§20 (derived
  quantities, application classes, the Cartesian boundary, time scaling). It
  is already referenced by three documents.
- **`ROADMAP.md`**: Stage 04 returns second derivatives; Stage 06 depends on
  Stage 04.
- **`README.md`**: compiler requirements if the minimum Clang changes; the
  `max_bandwidth` example after D5.
- **`docs/INDEX.md`**: a row for each new document.

## 8. Checks run

Python checks used numpy, scipy and sympy in a scratch virtual environment.
Compiler checks used the container's GCC 13.3, Clang 18.1.3, and Clang 19.1.1
with libc++ 18 (installed for this check). None of these files are in the
repository.

| Check | Result |
| --- | --- |
| `Bᵀ H B` bitwise symmetric, 10,000 random 3×3 | not symmetric in 9,330 |
| Pivoted-QR ratio vs κ₂, Kahan matrix 30×30 | 3.42 vs 9.88e3 |
| Pivoted-QR ratio vs κ₂, graded 40×10 | 9.764e8 vs 9.766e8 |
| Frobenius upper bound, same matrix | 9.89e8 |
| Singular values of `R` vs `A` | equal to rounding |
| 12 points on a circle, quadratic basis | rank 5 of 6 (random points: 6 of 6) |
| Collinear points, degree 1, before and after a random `B` | rank 2 both |
| `ω = −Δψ`, `W = −4 det H_ψ`, `Q = −W/4` (symbolic) | all hold |
| Support formulas near `q = 1` vs exact rationals | table in D8 |
| `sqrt(largest value below 1) < 1`, `float` and `double` | true for both |
| One-sided slope, linear vs quadratic fit | +0.05 vs ≈ 0 |
| Residual norm vs weighted RMS under weight scaling | table in D15 |
| NaN in one right-hand-side column | only that column non-finite |
| `std::expected` per toolchain | table in section 4 |
| tedder suite, Clang 19 and Clang 18 with libc++ (Debug) | 50/50 passed each |
