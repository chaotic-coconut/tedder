# tedder: field.hpp revision plan

Status: active

This plan implements the contracts recorded in
[ADR 0001](../../adr/0001-domain-periodic-contracts.md).

Keep the accepted scope: F1–F4, their regression tests, three small existing-type
tests, and one sanitizer job. Inspect the current checkout and retain subsequent
fixes to the reviewed
[c11b602 baseline](https://github.com/chaotic-coconut/tedder/tree/c11b60236c0d4634fefd5f72fe75b235e1f81fa6).
Keep the production changes in the existing header. Do not add kernels, fitting,
spatial indexing, I/O, dependencies, or new policy abstractions in this repair.

| Revision | Contract | Tests in PLAN.md |
| --- | --- | --- |
| F1 | Floating scalar, positive constant dimension, non-throwing operations; point representation stays flexible | B1 |
| F2 | Explicit periodic validity; invalid configurations admit no bandwidth | B2 |
| F3 | Guarded fast reduction and robust handling of finite outside-box coordinates | B3 |
| F4 | Accurate lifetime, indexing, numerical, and geometry documentation | A and documentation review |

Implement each behavioural revision with its matching tests. If a contract must
change, update both plans and the affected expectations together.

## F1. Tighten Domain without restricting point representation

Add `<type_traits>` and use a single `requires` block:

```cpp
template <class G>
concept Domain = requires(const G& g, const typename G::point_type& p,
                          const typename G::point_type& q) {
  typename G::scalar_type;
  typename G::point_type;
  // Check constant-expression suitability before evaluating dimension > 0.
  typename std::integral_constant<std::size_t, G::dimension>;
  requires Real<typename G::scalar_type>;
  requires G::dimension > 0;
  requires std::same_as<typename G::point_type,
                        Point<typename G::scalar_type, G::dimension>>;
  { G::dimension } -> std::convertible_to<std::size_t>;
  { g.offset(p, q) } noexcept -> std::same_as<typename G::point_type>;
  { g.distance_squared(p, q) } noexcept
      -> std::same_as<typename G::scalar_type>;
  { g.max_bandwidth() } noexcept -> std::same_as<typename G::scalar_type>;
};
```

The type requirement must precede `requires G::dimension > 0`; a separate helper
concept is unnecessary. Requires-expression checks stop in lexical order when
the result is determined. [C++ requires-expression rules](https://eel.is/c++draft/expr.prim.req)

Require `point_type` to equal `Point<scalar_type, dimension>`, per
[ADR 0002](../../adr/0002-point-type-restriction.md):

    requires std::same_as<typename G::point_type,
                          Point<typename G::scalar_type, G::dimension>>;

Place it after `requires Real<...>` and `requires G::dimension > 0`,
because `Point<scalar_type, dimension>` is only well-formed once both
hold. A model spelling `point_type` as `std::array<T,D>` still passes:
that is the same type.

The built-in methods already have `noexcept`; preserve them and the free helpers.
Reject potentially throwing domain methods at the concept boundary. Zero and
negative dimensions must fail cleanly, as must runtime dimensions. This excludes
zero-dimensional built-ins from `Domain` without banning `Point<T,0>` elsewhere.

## F2. Make periodic validity explicit

Keep aggregate construction and public mutable `length`. Add:

```cpp
constexpr bool is_valid() const noexcept;
```

Its result is true exactly when every period entry is finite. Finite positive
entries wrap; finite zero and negative entries remain open. Default construction
is valid and behaves as Euclidean geometry.

For valid configuration, `max_bandwidth()` remains half the shortest positive
period, or infinity when all axes are open. For invalid configuration it returns
zero, so `admits_bandwidth()` rejects every bandwidth. Keep strict positivity and
`h < max_bandwidth()`; equality is not allowed. Very small periods may leave no
representable positive admissible bandwidth.

Recompute validity and bandwidth from current `length`; do not cache them. Do not
add a throwing constructor, checked factory, or release-time checked geometry API.
`is_valid()` is specific to `Periodic`, not a new `Domain` requirement.

Document valid configuration and finite coordinates as preconditions of geometry
operations. Assert configuration validity at the start of `offset`; the existing
`distance_squared` path through `offset` then shares that debug check. Validation
queries themselves must work on invalid configurations without asserting.
Future reconstruction should validate configuration and bandwidth before its
neighbour loop; that API is outside this task.

## F3. Use the accepted finite, one-period fast path

Retain support for finite coordinates outside the canonical box. For each axis,
use the accepted structure below, with `L = length[d]`:

```cpp
scalar_type t = p[d] - q[d];
if (L > scalar_type{}) {
  if (std::abs(t) <= L) [[likely]] {
    t -= L * std::round(t / L);
  } else {
    t = std::remainder(p[d], L) - std::remainder(q[d], L);
    t -= L * std::round(t / L);
  }
}
r[d] = t;
```

The `std::isfinite(t)` test an earlier draft carried is unnecessary here.
A valid configuration has a finite period, so an infinite displacement
fails `std::abs(t) <= L` and routes to the slow path on its own. The
earlier draft needed it because it compared against `L * 2^(digits/2)`,
which could itself be infinite, making `inf <= inf` pass the guard.

This makes the guard depend on `is_valid()` holding. That is why `offset`
asserts it, and why the assertion is not redundant.

The fallback reduces each finite coordinate before subtracting. It therefore
handles the runtime case where the initial raw subtraction overflows. Open axes
retain ordinary subtraction and its representable-range limitation. Do not revive
the `L * 2^(digits/2)` threshold or route every normal in-box call through two
`remainder` calls.

Be precise about accuracy. In the fast branch the multiplier is only -1, 0, or
1; this avoids the large quotient/product cancellation defect. It does not make
the complete operation exact for arbitrary inputs: the initial `p-q` can already
have rounded. Neither branch recovers detail lost when coordinates were stored.
Do not promise universal bitwise equivalence between alternative formulations.

Preserve these semantics:

- `offset(p,q)` means `p-q`, reduced to a minimum image on periodic axes.
- For period 2, `offset({0},{1}) == {1}` and the reverse is `{-1}`.
- Retain antisymmetry within floating-point tolerance. At other half-period
  ties, the sign follows the selected path's reduced displacement and
  `std::round` convention. Equivalent coordinate representations can select
  opposite equally valid images; do not promise invariant signs at those ties.
- Bandwidth remains strictly below the half-period limit. This does not remove
  the public geometry operation's existing canonical tie contract.
- The four independently specified regressions in PLAN.md B3 must pass.

`[[likely]]` may remain, but it is a compiler hint, not a guarantee of code layout
or speed. The supplied benchmark motivates this design; its percentages are not
a portable acceptance criterion. [C++ likelihood attribute](https://eel.is/c++draft/dcl.attr.likelihood)

Floating-point exception flags/trapping modes are not a new guarantee of this
repair. In particular, recovering a finite result after raw subtraction overflow
does not promise that no floating-point exception flag was raised.

## F4. Correct documentation without expanding the API

### Storage and access

- Fix `SampleView`'s lifetime wording: the underlying point/value storage must
  stay alive and valid when accessed. The original span descriptors need not
  outlive the view; they are copied. Owner destruction or reallocation can
  invalidate the view.
- Keep equal lengths as a constructor precondition with a debug assertion.
  Ordinary indexing remains unchecked. Document `i < size()`, `r < rows`, and
  `c < cols`; for subviews use `first <= size()` and
  `count <= size() - first`.
- Preserve `LocalFit<T,D,C>::jacobian` as `Matrix<T,C,D>`, where `[c,d]` means
  the derivative of component `c` along coordinate `d`.
- Describe `is_flat` as a storage-size observation. It alone does not establish
  portable serialization or permission to treat nested arrays as one scalar
  array. Do not add I/O code or lifetime restrictions to satisfy deferred tests.

### Geometry and numerical range

Attach isotropy wording to the built-in Euclidean/Periodic geometries, not to
`Domain`: custom models need not use the same metric. A bandwidth spanning a
thin direction reduces locality there; it does not by itself remove a dimension.
Full rank still depends on the sample layout and polynomial degree, not merely
on coordinates varying along every axis.

Treat axis rescaling as an optional modelling choice because it changes the
metric and neighbours. If `u_d = (x_d-a_d)/s_d`, convert derivatives with
`J_x[:,d] = J_u[:,d]/s_d`. Scaling or centring values does not improve the design
matrix's conditioning, but can help manage numerical range. Do not prescribe
normalisation by domain extent as a universal preparation step.

Explain that keeping periodic coordinates in a canonical box normally uses the
fast path; it is a recommendation, not an input restriction or an exactness
guarantee. Put extended data-preparation guidance in `PROJECT.md`. Keep README
short and limited to current capabilities and the next few features.

Squared distances can overflow or underflow even when offsets are finite.
`distance()` currently takes `sqrt(distance_squared())`, so it can lose a
representable final norm through an intermediate result. Document this limitation;
do not add a stable-norm API or change generic distance to assume an isotropic
metric. Avoid fixed overflow thresholds that ignore dimension and coordinate
magnitudes.

### Constant evaluation

`std::round` and `std::remainder` are constexpr in C++23; `std::sqrt` follows in
C++26. Compiler and standard-library support must be distinguished from the
standard requirement. An absent feature macro alone does not prove that every
relevant operation is unavailable. [P0533R9](https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p0533r9.pdf),
[P1383R2](https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2023/p1383r2.pdf)

Keep valid constexpr declarations and existing supported static checks. Do not
require periodic offset evaluation at compile time on every current toolchain,
or add a capability-check framework in this repair. Runtime tests remain required;
runtime overflow recovery is not a promise of constant evaluation for that input.

## Acceptance and handoff

F1 fixtures reject unsupported domains cleanly while accepting a custom point
representation. F2 rejects non-finite periods and reflects mutations. F3 passes
the four regressions and preserves existing geometry signs, open axes, and strict
bandwidth limits. The three A tests pass without redefining their existing APIs.

Run the Debug/Release and sanitizer checks in PLAN.md, reporting unavailable
configurations honestly. Summarise API/contract changes and remaining numeric
limits. The consumer project, assertion subprocesses, hard compile-fail harness,
two-translation-unit check, exhaustive type matrix, and `long double` coverage
remain deferred as agreed.
