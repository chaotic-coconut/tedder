# ADR 0004 — local fit results, diagnostics and errors

Status: proposed. Draft for the author's review; not authoritative until
accepted. Depends on [ADR 0003](0003-geometry-local-frame-support.md).

Implementation: the `LocalFit` result type (§1 and §2) is already on `main`,
in commit e99f6a5. Everything else is for Stage 03 (degree 0 and 1 fitting)
and later.

## Context

Before e99f6a5, `LocalFit` held a value, a Jacobian, `neighbors` and
`bandwidth`. No fitting code exists yet.

The implementation plan of 2026-10-04 and its
[review](../reviews/2026-10-04-implementation-plan-review.md) agreed on the
following:

- separate the polynomial degree from the derivatives the result holds;
- return derivatives only in physical coordinates;
- keep diagnostics outside the result;
- keep solver evidence separate from the acceptance policy;
- use Eigen for linear algebra;
- use `std::expected` for errors.

The author then reversed two of the review's recommendations. Non-finite
sample values are rejected instead of dropped, and the return type is
`std::expected<FitSuccess, FitFailure>` instead of a separate outcome struct.
The author also chose the name `FitOrder` and the default order `linear`,
both implemented in e99f6a5.

## Decision

### 1. Result order and polynomial degree

- **Result order.** `FitOrder` says which derivatives a result holds:

  ```cpp
  enum class FitOrder : unsigned { constant = 0, linear = 1, quadratic = 2 };

  template <Real T, std::size_t D, std::size_t C, FitOrder P = FitOrder::linear>
    requires(P <= FitOrder::quadratic)
  struct LocalFit;
  ```

  `constant` holds the value; `linear` adds the Jacobian; `quadratic` adds the
  Hessians. Because the constraint is on the primary template, naming an
  order above `quadratic` fails where the type is named. Without the
  constraint, the result would only be an incomplete type, which can still be
  named.
- **Default `linear`.** The cheapest fit, degree 0, is a weighted average and
  does not reproduce linear fields. With `f = 2 + 3x` and `h = 0.1`, it erred
  by 9e-3 to 5.5e-2 wherever the neighbourhood was not symmetric (uneven
  sampling, the edge of the data). Degree 1 was exact to rounding. Once the
  polynomial has degree 1, returning its Jacobian costs only `C·D` divisions.
- **Degree is separate.** Fitting functions take the polynomial degree `n`
  separately, constrained so that `n` is at least the order. At first only
  `n` equal to the order is accepted. Allowing one degree higher later must
  not change the result type. One example is a quadratic fit returning a
  `linear` result, for a better gradient at the edge of the data.
- **Naming caveat.** The enumerators read like polynomial degrees, but they
  describe what the result holds. If a higher degree is allowed later, a
  `FitOrder::linear` result could come from a quadratic fit. Revisit the
  names at that point (see open points).
- **Generic code.** `WithJacobian<P>` and `WithHessians<P>` are named concepts,
  with `WithHessians` built from `WithJacobian`.
  - Overloads constrained by them are ordered by subsumption, so a quadratic
    result picks the Hessian overload.
  - Two overloads constrained by the raw comparisons `(P >= linear)` and
    `(P >= quadratic)` are ambiguous for a quadratic result, on GCC 13 and on
    Clang 19.
  - The specializations do not inherit from one another; generic code
    constrains on the order instead.

### 2. Result contents

`LocalFit` holds the reconstructed field at `x₀`, in physical coordinates
only.

| `FitOrder` | Members |
| --- | --- |
| `constant` | `value`: `Value<T, C>` |
| `linear` | adds `jacobian`: `Matrix<T, C, D>`, with `jacobian[c, d] = ∂F_c/∂x_d` (unchanged) |
| `quadratic` | adds `hessians`: `std::array<Matrix<T, D, D>, C>`, with `hessians[c][i, j] = ∂²F_c/∂x_i∂x_j` |

- **Removed fields.** `neighbors` and `bandwidth` are removed. Counts are
  diagnostics (§6), and the bandwidth is an input the caller already has.
- **Nothing else.** The result holds no residuals, ranks, condition estimates,
  counts or derived quantities such as divergence or vorticity.
- **Bitwise symmetry.** Every fit produces bitwise symmetric Hessians:
  `hessians[c][i, j] == hessians[c][j, i]`.
  - This is a duty of the fitting code, not of the type. A `LocalFit` filled by
    hand can be asymmetric.
  - The physical transform `H_x = Bᵀ H_ξ B` (ADR 0003) computes the entries
    with `i ≤ j` and copies them to `j > i`.
  - Computed independently, the two triangles can differ by a few ulps; this
    happened in 9,330 of 10,000 random 3×3 cases.
  - For `B = I/h`, symmetry holds without the copy.

### 3. Basis

- **Space.** A complete total-degree basis in the local coordinates `ξ` of
  ADR 0003.
- **Taylor normalisation.** `p(ξ) = Σ_{|α| ≤ n} c_α ξ^α / α!`, so each
  coefficient `c_α` equals `∂^α p(0)`, the ξ-derivative at `x₀`. For `n ≤ 1`
  every `α! = 1`.
- **Ordering.** Terms are ordered by total degree. Degree 1 terms follow axis
  order. Degree 2 terms are pairs `(i, j)` with `i ≤ j`, in lexicographic
  order. For `D = 2, n = 2` that gives:

  ```text
  1, ξ₀, ξ₁, ξ₀²/2, ξ₀ξ₁, ξ₁²/2
  ```

- **Basis size.** `m = C(D+n, n)`: 1, `D + 1`, `(D + 1)(D + 2) / 2`. It is
  exposed as `constexpr`, not stored.
- **Physical derivatives.** Derivatives reach physical coordinates only
  through the frame transform of ADR 0003. Derived quantities never use
  ξ-derivatives.

### 4. Linear system and solver

- **Design matrix.** Each fit builds one weighted design matrix `A = W^{1/2} Φ`.
  - Rows are the active samples, in ascending sample index (ADR 0003 §5).
  - Columns are the basis terms.
  - `A` depends on geometry, the frame and the kernel weights, never on the
    values.
- **One factorisation.** All `C` right-hand sides are solved against one
  factorisation of `A`.
- **Seams for future observations.** These keep derivative observations
  possible later, without promising that extraction or result types would stay
  unchanged.
  - Row construction is separate from factorisation.
  - The solver does not assume one row per sample.
  - It accepts rank-deficient systems and reports their rank.
  - `LocalFit` is assembled from the coefficients in a separate step.
- **Solver roles.** Column-pivoted QR is the default. SVD serves as the test
  reference and as an optional diagnostic. A normal-equations path is added
  only after the benchmarks in plan §22. There is no automatic fallback.
- **Eigen.** Eigen is the linear-algebra library, with these limits:
  - Eigen appears only in namespace `detail`, never in a public signature.
  - `field.hpp` stays free of Eigen; fitting lives in its own header.
  - tedder never defines `EIGEN_USE_*` macros. Routing to BLAS, LAPACK or MKL
    is the consumer's choice.
  - The design matrix is stored column-major internally.
  - There is no backend abstraction until a second backend exists.
- **Numerical rank.** Using the R factor of the column-pivoted QR, the rank is
  the number of diagonal entries with `|r_kk| > τ |r_11|`. `τ` is a documented
  multiple of `ε(T)`, chosen from the degeneracy experiments. It is tedder's
  own definition, not a library default.

### 5. Solver evidence and conditioning

The solver reports facts. The policy in §7 decides what they mean.

- **Rank.** The numerical rank as defined in §4.
- **Condition estimate.** Optional, and always labelled with its provenance:

  | Provenance | Value | Meaning |
  | --- | --- | --- |
  | `cpqr` | `|r_11 / r_mm|` | A lower bound on `κ₂(A)` that can be far too low (Kahan 30×30: 3.42 against 9.88e3). Never reported as `κ₂`. |
  | `svd` | `σ_max / σ_min` of the `m × m` factor `R` | An estimate of `κ₂(A)`. `A` and `R` have the same singular values, and `m ≤ 10`. Accuracy falls as `κ₂` approaches `1/ε(T)`. |
  | none | — | No estimate computed. |

  Measured in `double` on matrices built with a given `κ₂`, the `svd`
  estimate's relative error was `5e-14` at `κ₂ = 1e4` and `4e-2` at `1e15`.
- **Full SVD of `A`.** A test reference only, not used in production.
- **Normal equations.** If ever added, their estimates get their own
  provenance, because they describe `AᵀA` rather than `A`.

### 6. Diagnostics

Diagnostics are kept outside `LocalFit`. They contain:

- **Counts.** `candidates` (after duplicate removal and the support test of
  ADR 0003), and `active` (candidates with a finite, positive weight).
- **Solver evidence** from §5.
- **Residual quality, per component.**

  ```text
  weighted_residual_rms = sqrt(Σ wᵢ eᵢ² / Σ wᵢ),   eᵢ = p(ξᵢ) − yᵢ
  ```

  The sums run over the active samples. The value does not change when all
  weights are rescaled together. It measures fit quality. Kernel weights
  localise the fit and are not inverse variances, so it is not a calibrated
  estimate of measurement noise.
- **Degrees of freedom.** `active − numerical rank`, derived rather than
  stored. At zero, the residual is zero by construction and carries no
  information.

There are two diagnostic types, which share their sub-structures (counts,
solver evidence, residuals):

- **Success diagnostics.** Every member is present; nothing is optional.
- **Failure diagnostics.** One optional member per stage, present exactly when
  that stage ran.

### 7. Errors, policy and return type

```cpp
// Conceptual; names and layout are open.
std::expected<FitSuccess, FitFailure>

struct FitSuccess { LocalFit<T, D, C, P> fit;  /* success diagnostics */ };
struct FitFailure { FitError error;            /* failure diagnostics */ };
```

Policy for value fitting in Stages 03 and 04. The checks run in table order,
and the first one that fails is reported.

| `FitError` | Stage | Condition | The failure carries |
| --- | --- | --- | --- |
| `inadmissible_bandwidth` | Validation | `!admits_bandwidth(g, h)`. This includes an invalid `Periodic`, whose `max_bandwidth()` is 0. | Nothing |
| `nonfinite_evaluation_point` | Validation | A coordinate of `x₀` is not finite | Nothing |
| `nonfinite_sample_value` | Filter | An active sample has a non-finite value component | Counts, the number of such samples, and the index of the first one |
| `insufficient_samples` | Filter | `active < m` | Counts |
| `rank_deficient` | Factorisation | Numerical rank `< m` | Counts, solver evidence |
| `nonfinite_result` | Solve | A coefficient or derivative is not finite | Counts, solver evidence |

Policy rules:

- **Rank.** Success requires full numerical rank `m`. This rule belongs to
  the fitting policy, not the solver; a future observation model may accept
  rank deficiency under a different policy.
- **Conditioning.** A high condition number is not an error. Callers may
  reject fits using the diagnostics.
- **Non-finite sample values.**
  - Only active samples are checked. Values at zero-weight or out-of-support
    samples do not matter.
  - They are rejected by default, and no sample is ever dropped silently.
  - A drop policy may be added later as an explicit option, with rejection
    remaining the default. It would report dropped samples in the counts.
- **Non-finite sample coordinates.** They remain a precondition (ADR 0001)
  and are not detected; see ADR 0003, open points.

Recorded statement on rank:

> For a fixed active sample set with strictly positive weights and a complete
> total-degree basis, an invertible affine change of the polynomial
> coordinates preserves the exact rank of the design matrix. Conditioning, and
> therefore the numerical rank, can change. Because `B` also determines
> support and weights, changing `B` can change the active set and therefore
> the exact rank.

Toolchain:

- **Requirement.** `std::expected` requires `__cpp_lib_expected`.
- **Tested on 2026-10-04 (Ubuntu 24.04 packages):**
  - GCC 13.3 provides it.
  - Clang 18.1.3 with libstdc++ 13 does not. libstdc++ requires
    `__cpp_concepts ≥ 202002L`, and Clang 18 defines `201907L`.
  - Clang 18 with libc++ 18 provides it, and the existing suite passes.
  - Clang 19.1.1 with libstdc++ 13 provides it, and the existing suite passes.
- **CI.** Clang 18 is replaced by Clang 19, in both the build matrix and the
  sanitizer job.
- **Readable failure.** The fitting header includes `<version>`, checks
  `__cpp_lib_expected`, and stops with a readable `#error`.
- **README.** It states that the fitting header is unavailable with Clang 18
  and libstdc++. There is no shim and no formal minimum compiler version.
- **`field.hpp`.** It does not use `std::expected`.

### 8. Floating-point types

- **Supported types.** `float` and `double` are supported and tested for
  every fitting feature.
- **Tolerances.** Every tolerance and rank threshold is a multiple of
  `std::numeric_limits<T>::epsilon()`. No absolute literals.
- **Degeneracy sweeps.** They run down to about `ε(T)`, which is 1e-1 to
  1e-7 for `float`.
- **Test accuracy.** Reproduction-test tolerances scale with the conditioning
  of the problem. When `κ₂` is large, `float` Hessians may keep only a few
  correct digits.
- **Compute type.** An internal `compute_type` alias equals `T`. Changing it,
  for example to `double` for `float` input, is a later decision.
- **Range.** Scaling by `1/h²` overflows `float` when `h` is below about
  `5e-20`. This is documented, not guarded.

## Consequences

Done for §1 and §2 in commit e99f6a5:

- **`field.hpp`.** `FitOrder`, the constrained `LocalFit` primary template with
  three specializations, `WithJacobian` and `WithHessians`. `neighbors` and
  `bandwidth` are removed.
- **`tests/test_types.cpp`.** These tests cover the result type, for `double`
  and `float`:
  - the default order is `linear`;
  - each order has the right members;
  - the concepts give the right truth table;
  - an order above `quadratic` is rejected, checked through a concept rather
    than a compile-fail harness;
  - Hessian shapes with `D ≠ C`, every entry independently addressable, and
    value-initialised zeros.

Still to do:

- **New fitting header.** It depends on Eigen; how CMake finds Eigen is an
  open point.
- **CI and README.** The CI and README changes from the toolchain list in §7.
- **Stage 03 and 04 files.** They carry the acceptance tests, including:
  - constant, affine and quadratic reproduction, for scalar fields and for
    `C ≠ D`;
  - Hessian symmetry checked with exact `==`, including one non-diagonal
    internal frame;
  - degenerate layouts that fail with `rank_deficient`:
    - collinear points (2D, degree 1);
    - coplanar points (3D);
    - points on any conic for degree 2: an off-centre circle, a rotated
      ellipse, two parallel lines, and a plus-shaped stencil (each rank 5 of
      6);
    - points on a sphere in 3D (rank 9 of 10);
  - a non-finite value at an active sample gives `nonfinite_sample_value`; at a
    zero-weight or out-of-support sample the fit succeeds;
  - `float` and `double` throughout.
- **`ROADMAP.md`.** Stage 04 returns Hessians. Stage 06 depends on Stage 04,
  because the Laplacian needs Hessians.
- **On acceptance.** Mark this ADR accepted in `docs/INDEX.md`, and link it
  from the Stage 03 and Stage 04 rows of `ROADMAP.md`.

## Rejected alternatives

- **Defining `LocalFit`'s parameter as the polynomial degree.** It would rule
  out a higher-degree fit returning fewer derivatives.
- **Default order `constant`.** Cheapest, but biased wherever the
  neighbourhood is not symmetric (§1).
- **Inheritance between the `LocalFit` specializations.** It would allow
  slicing and complicate aggregate initialisation. Generic code constrains on
  the order instead.
- **Diagnostics inside `LocalFit`.** That mixes reconstructed data with solver
  facts.
- **An unconditional `FitFailure::numerical_rank`.** No rank exists before
  factorisation.
- **A separate outcome struct, `{ expected<LocalFit, FitError>, diagnostics }`.**
  Success diagnostics would need optional members, so "complete on success"
  would be a documented promise rather than a property of the type.
- **Dropping non-finite samples by default.** A caller who ignores the
  diagnostics would receive a successful fit of different data.
- **Reporting the pivoted-QR diagonal ratio as `κ₂`.**
- **Describing the SVD value as exact.**
- **A compatibility shim for `std::expected`.**

## Open points

- **Names.** `FitSuccess`, `FitFailure` and the diagnostic structs. Also
  whether the `FitOrder` enumerators should describe derivatives rather than
  degrees, which becomes pressing once a higher degree than the order is
  allowed (§1).
- **`τ`.** Its value for each scalar type, from the degeneracy experiments.
- **Eigen in CMake.** The minimum Eigen version, and whether to use
  `find_package(Eigen3 CONFIG)`, `FetchContent`, or both.
- **`compute_type` for `float`.**
- **Optional diagnostics.**
  - Whether to add the effective sample size `(Σw)² / Σw²`.
  - Whether the `svd` condition estimate is computed by default or on request.
- **Drop policy.** An explicit drop policy for non-finite values, if a real
  application needs one.
- **Out of scope.** Derived quantities (Stage 06) and the content of
  `PROJECT.md`.
