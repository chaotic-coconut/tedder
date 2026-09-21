# tedder — decisions on the field.hpp revision (F1–F4)

Status: accepted, with corrections. The `point_type` decision is superseded
by [ADR 0002](0002-point-type-restriction.md). The F3 fast-path guard was
later simplified; see
[field-plan.md](../stages/00-foundation/field-plan.md) F3 for the shipped
form. Two claims in the body are corrected in the Corrections section at
the end. All other decisions stand.

Measurements on GCC 13.3 and Clang 18, x86-64 Linux.

## Summary

| ID | Decision |
| --- | --- |
| F1 | Accept, minus the `point_type` restriction. Single `requires` block, type requirement first. |
| F2 | Accept as proposed. |
| F3 | Accept the defect. Fast path restricted to the provably exact regime, plus a finiteness test. |
| F4 | Accept, with the isotropy claim corrected and the rescaling advice softened. |

Three small tests recovered from the deferred set. The rest of the companion
test plan stays deferred.

### Rejected approaches

1. The `2^(digits/2)` threshold was never exact. Counterexample below.
2. The threshold had an overflow hole that returned NaN.
3. A separate `ConstantDimension` concept is not required; lexical ordering
   inside one `requires` block is sufficient and is guaranteed by the standard.
4. `Euclidean` and `Periodic` already have `noexcept` on all three methods. No
   change needed there.
5. A bandwidth spanning a thin axis does not make the reconstruction
   two-dimensional. It removes locality along that axis. The fit stays full rank.

## Confirmed defects

All reproduced against the baseline header.

```
F1  Domain<StringScalar> = 1     scalar_type = std::string satisfies Domain
F1  Domain<MutableDim>   = 1     a runtime `static int dimension` satisfies it
F1  Domain<Throwing>     = 1     and the noexcept free helpers then terminate

F2  infinite period: max_bandwidth = inf, admits_bandwidth(1.0) = true,
                     but distance_squared(a, a) = NaN
F2  NaN period:      silently behaves as an open axis

F3  double L=3, offset(1e16, 0) = -2      correct value is 1
F3  float  L=3, offset(1e8,  0) =  0      correct value is 1
```

F3 corrects an earlier claim that the formula was sound and the information was
lost in the input. That was wrong: `1e16` is exactly representable and
`std::remainder(1e16, 3)` returns exactly 1. The loss happens in
`L * round(t/L)`, where the product rounds at that magnitude.

## F1. Tighten Domain

### Accepted

- `scalar_type` must satisfy `Real`.
- `dimension` must be a positive compile-time constant.
- `offset`, `distance_squared` and `max_bandwidth` must be `noexcept`, since the
  free helpers that call them are unconditionally `noexcept` and would otherwise
  terminate.

The built-in models already satisfy all three. This closes the door on future
models, it does not change existing code.

### Rejected

`point_type` must be exactly `Point<scalar_type, dimension>`. That removes the
only freedom the concept has, for nothing the project needs today. A design
narrowing, not a repair. Revisit when a concrete need appears.

### The constant-dimension check

A plain nested requirement is not enough. A non-constant `dimension` produces a
hard error rather than a `false` result:

```
error: the value of 'MutableDim::dimension' is not usable in a constant expression
```

Nested requirements do not SFINAE on "not usable in a constant expression". The
value has to appear in a deduced context, which a template argument provides.

`[expr.prim.req]` states that substitution and semantic constraint checking
proceed in lexical order and stop at the first unsatisfied condition, so the type
requirement placed first inside the same block is sufficient. No separate concept
is needed.

```cpp
template <class G>
concept Domain = requires(const G &g, const typename G::point_type &p,
                          const typename G::point_type &q) {
  typename G::scalar_type;
  typename G::point_type;
  // Forces dimension into a template argument, so a non-constant one is a
  // substitution failure rather than a hard error. Must come first.
  typename std::integral_constant<std::size_t, G::dimension>;
  requires Real<typename G::scalar_type>;
  requires G::dimension > 0;
  { G::dimension } -> std::convertible_to<std::size_t>;
  { g.offset(p, q) }           noexcept -> std::same_as<typename G::point_type>;
  { g.distance_squared(p, q) } noexcept -> std::same_as<typename G::scalar_type>;
  { g.max_bandwidth() }        noexcept -> std::same_as<typename G::scalar_type>;
};
```

Rejects a negative dimension for free: `integral_constant<size_t, -1>` is a
narrowing conversion in a template argument, so substitution fails.

Verified on both compilers, every rejection clean and no hard errors:

```
Good 1    StringScalar 0    MutableDim 0    ZeroDim 0    NegDim 0
Throwing 0    int 0
```

Requires `<type_traits>`.

Side effect: `Euclidean<T,0>` and `Periodic<T,0>` no longer model `Domain`. The
`Point<T,0>` alias is unaffected.

## F2. Invalid periodic configuration

Accepted as proposed, unchanged.

- `constexpr bool is_valid() const noexcept`, true exactly when every entry in
  `length` is finite.
- Invalid configuration makes `max_bandwidth()` return zero, so the existing
  strict positive-and-below-limit rule rejects every bandwidth.
- Finite zero or negative periods continue to mean open axes. Default
  construction stays valid and equivalent to Euclidean.
- `length` stays public and mutable, so validity is computed from current
  contents and never cached.
- `offset` and `distance_squared` require a valid configuration: document the
  precondition, add a debug assertion, no release-time checked access.

No exceptions, no throwing constructor, no owning domain type. `is_valid()` is
specific to `Periodic` and is not required of every `Domain`.

## F3. Robust periodic reduction

### Why the first proposal failed

A rejected fast path guarded by
`|t| <= L * 2^(digits/2)` and claimed exactness below that limit. Both the claim
and the guard were wrong.

**Accuracy.** With `L = nextafterf(1.0f, 2.0f)`, `p = 4000.25f`, `q = 0.0f`, all
exactly representable, the fast path is selected and returns
`0.24951171875` where the correct reduction is `0.249523162841796875`. Error
1.14e-5, which is 1.14e-5 of a period. The threshold was a rule of thumb whose
error bound was never computed; it admits errors on the order of
`L * 2^(-digits/2)`.

**Overflow.** With `p = 1e308`, `q = -1e308`, `L = 1e308`, both `|p - q|` and
`L * 2^(digits/2)` evaluate to infinity. `inf <= inf` is true, the fast path is
taken, and `offset` returns NaN where the correct answer is 0.

### Why unconditional `remainder` is also rejected

Measured over 4.2M in-box pairs, `double`, `L = 3`:

| formulation | time |
| --- | ---: |
| original, unguarded | 31.2 ms |
| guarded, `|t| <= L` and finite | 38.6 ms |
| always `remainder` (two calls per axis) | 199.4 ms |

Six times slower in the innermost loop of the library. The guard costs 24%.

### The accepted guard

Drop the magic constant. Use the regime that is provably exact: when the
displacement is finite and already within one period, `round(t/L)` is `-1`, `0`
or `1`, the product `L * k` is exact, and Sterbenz's lemma makes the subtraction
exact. Every other case, including non-finite, reduces each coordinate first.

```cpp
scalar_type t = p[d] - q[d];
if (length[d] > scalar_type{})
  {
    // Exact when the displacement is already within one period. Anything else,
    // including a non-finite difference, reduces each coordinate first.
    if (std::isfinite(t) && std::abs(t) <= length[d]) [[likely]]
      t -= length[d] * std::round(t / length[d]);
    else
      {
        t = std::remainder(p[d], length[d]) - std::remainder(q[d], length[d]);
        t -= length[d] * std::round(t / length[d]);
      }
  }
r[d] = t;
```

Verification:

```
counterexample 1  0.249523162841796875, matches remainder bit for bit
counterexample 2  0                                   (was NaN)
regression 1e16   1                                   (was -2)
regression 1e8f   1                                   (was 0)
2M random pairs   worst |error| / L = 5.6e-17, one ulp
```

The 5.6e-17 residual appears only in the slow path, where both formulations
perform the same final reduction. The fast path is exact in its regime, with no
qualification.

`[[likely]]` is kept with the understanding that it buys nothing at runtime on a
perfectly predicted branch. It moves the cold path out of line, which helps
instruction cache.

If the 24% ever matters, the honest lever is an explicit in-box precondition on
the domain, not a looser threshold.

### Preserved

- `offset(p, q)` remains the displacement `p - q` reduced to the minimum image.
- The canonical tie for `L = 2`: `offset({0},{1}) == {1}`, reverse `{-1}`.
- Strict exclusion of the half-period boundary in `admits_bandwidth`.
- Support for finite coordinates outside the canonical box.

### Not fixed by F3

`distance_squared` still squares the reduced offset, so a domain with enormous
periods can overflow the square even when `offset` is correct. Documented as a
numeric-range limitation rather than repaired. The sum of squares overflows above
roughly 7.7e153 per coordinate in `double` and 1.1e19 in `float`.

## F4. Comments and documentation

The accepted guidance corrects the isotropy claim and treats coordinate
rescaling as a modelling choice.

### Isotropy, corrected

An earlier analysis incorrectly assumed a bandwidth spanning a thin axis makes the
reconstruction "silently two-dimensional". That is wrong. Sample coordinates
along the thin axis still vary within the neighbourhood, so the design matrix
stays full rank and the derivative along that axis is still estimated. What is
lost is *locality* along that axis: the fit degenerates toward a global
polynomial in that direction.

For a domain 1000 units wide and 2 tall:

| h | fraction of x covered | fraction of z covered |
| ---: | ---: | ---: |
| 1 | 0.1% | 50% |
| 10 | 1.0% | 100% |
| 100 | 10% | 100% |

On `Domain`:

```cpp
// Distances are isotropic: the same bandwidth applies to every axis. With very
// different extents per axis, a bandwidth chosen for the wide direction spans
// the whole thin one, so there is no locality left in it.
```

On `Periodic`:

```cpp
// Coordinates should already be inside the box. Points far outside still work,
// but the reduction is exact only for a displacement within one period.
```

### Rescaling, softened

Rescaling axes is a modelling choice, not a correction. It changes the metric
and therefore which points are neighbours. Aspect ratio alone does not justify
it.

README section:

```markdown
## Preparing your data

`tedder` treats all axes the same way. One bandwidth applies in every direction,
so the shape of your domain affects what counts as a neighbourhood.

**Axes with very different extents.** A bandwidth suited to a wide direction may
span a thin one entirely, leaving no locality along it. Rescaling the axes is one
way to handle this, but it changes the metric and therefore the neighbourhoods,
so treat it as a modelling decision rather than a fix.

If you rescale with `u_d = (x_d - a_d) / s_d`, convert the gradient back with

    J_x[:, d] = J_u[:, d] / s_d

**Wrap periodic coordinates into the box.** Points outside still work, but the
reduction is exact only for a displacement within one period.

**Values.** Scaling or centring values does not improve the conditioning of the
design matrix, since it does not appear there. It can still help keep the
numerical range manageable.
```

The gradient conversion is the same class of factor as `J_j = b_j / h` in the
project brief, and it goes wrong the same way if forgotten.

### constexpr correction

`std::round` is constexpr in C++23 (P0533); `std::sqrt` only in C++26 (P1383). My
earlier claim that both wait for C++26 was wrong about the standard.

The practical conclusion survives on different grounds: neither GCC 13 nor Clang
18 defines `__cpp_lib_constexpr_cmath`, so neither standard library implements
P0533. GCC accepts `constexpr std::round` through its own builtin folding; Clang
rejects it. Use a capability check if compile-time periodic evaluation is wanted,
and keep runtime coverage regardless.

## Test scope

### Accepted now

Phase B regression tests tied to F1, F2 and F3, including both of
counterexamples as named cases.

Three tests recovered from the deferred set. Each exercises tedder's own indexing
and wiring rather than the standard library, and each is cheap:

- rectangular `Matrix`, checking row-major indexing
- `LocalFit` with `C != D`, checking the Jacobian's extents
- `SampleView::subview`, checking points and values stay aligned

The ASan/UBSan CI job.

### Still deferred

Most of the exhaustive `Matrix` and `SampleView` coverage, the child-process
assertion test, compile-fail checks, the separate consumer CMake project, the
two-translation-unit link test, and `long double` coverage. These belong to a
mature library. `tedder` cannot yet reconstruct a field.

Checking that mutating an owner is visible through a `SampleView` tests
`std::span`, not this code.

## Corrections

**Fast-path exactness.** The body says: "The fast path is exact in its regime,
with no qualification." That holds for the wrapping step, not for the whole
`offset` operation, because the initial subtraction `p - q` can round before
the reduction runs. Counterexample, binary64: `L = 1`, `p = 0x1p-54`, `q = 1`.
The difference rounds to `-1`, the fast path returns `0`, and the exact
minimum image is `0x1p-54`. The correct claim: for an already-rounded finite
displacement `t` with `|t| <= L`, the wrapping correction introduces no
further rounding. When `round(t/L)` is zero `t` is unchanged; when it is ±1
the product `L * k` is exact and the subtraction falls under Sterbenz's
lemma.

**Full rank.** The body says: "Sample coordinates along the thin axis still
vary within the neighbourhood, so the design matrix stays full rank and the
derivative along that axis is still estimated." The same claim appears in the
summary list near the top: "The fit stays full rank." The correction applies
to both. Variation per axis is not
sufficient: collinear points `(0,0), (1,1), (2,2)` vary in both coordinates,
yet the degree-one design matrix `[1, x, y]` has rank 2. What a bandwidth
spanning a thin axis loses is locality along it. Rank deficiency is a
separate question the fit must detect.
