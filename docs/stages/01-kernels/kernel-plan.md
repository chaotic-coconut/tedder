# tedder: Stage 01 kernel plan

Status: draft for the author's review. It becomes the active Stage 01 plan
once [ADR 0003](../../adr/0003-geometry-local-frame-support.md) is accepted.

This plan implements the kernel contract of ADR 0003 §4: the kernel takes
`q`, its precondition is `0 ≤ q ≤ 1`, `w(1) = 0`, the result is finite and
non-negative, it is positive inside the support in exact arithmetic, it
carries no normalising constant, and it is equality-preserving and free of
side effects. The ADR is authoritative; this plan does not repeat its
reasons.

## Scope

In scope:

- `include/tedder/kernel.hpp`: the `Kernel` concept and the Stage 01 kernels.
- `tests/test_kernel.cpp`, registered in `tests/CMakeLists.txt`.

Out of scope:

- the local frame, `q²` and the candidate test (Stage 02);
- neighbour search (Stage 02);
- weights inside a fit, and counting active samples (Stage 03);
- kernel derivatives, and kernels that are not compactly supported;
- kernels other than Wendland C2 and C4. Wendland C6 and Epanechnikov are
  possible later additions, and neither would change the concept.
  Epanechnikov, `1 − q²`, would need the factored form `(1 − q)(1 + q)`.
  Computed as `1 − q·q`, its relative error near `q = 1` reached 724
  `epsilon` in `float` and about 1.1e7 `epsilon` in `double` in a check
  against exact values; the factored form stayed below 1 `epsilon`.

`field.hpp` does not change.

## Order of work

The workflow guide asks for a concrete kernel before the concept, and for a
second model before the concept settles:

| Step | Work | Tests |
| --- | --- | --- |
| K1 | Wendland C2 as a concrete kernel | T2 for C2 |
| K2 | The `Kernel` concept, checked against C2 | T1 |
| K3 | Wendland C4, the second model | T1 and T2 for C4 |

K3 tests the concept: C4 must satisfy it unchanged. If it needs a change to
the concept, stop and record why before changing it.

The author designs and writes K1 to K3. AI may generate the T1 and T2 tests
afterwards, when delegated.

## K1 and K3. The kernels

Wendland's functions for dimension 3 with smoothness index 1 and 2, scaled
so that `w(0) = 1` (H. Wendland, Advances in Computational Mathematics 4,
1995, 389–396):

| Kernel | `w(q)` on `[0, 1]` | Smoothness of `x ↦ w(‖x‖)` | Near `q = 1` |
| --- | --- | --- | --- |
| C2 | `(1 − q)⁴ (4q + 1)` | C² | `w(1 − δ) ≈ 5 δ⁴` |
| C4 | `(1 − q)⁶ (35q² + 18q + 3) / 3` | C⁴ | `w(1 − δ) ≈ (56/3) δ⁶` |

Wendland's functions are positive definite on `Rᵈ` for `d ≤ 3`. That
matters for RBF interpolation, not here. As least-squares weights they serve
any `D`.

Checked symbolically with SymPy, for both kernels:

- `w(0) = 1` and `w(1) = 0`;
- no root in `[0, 1)`, so `w > 0` there (ADR 0003 §4, item 4);
- the expansion at 0 has no odd powers below `q³` (C2) or `q⁵` (C4), which
  is what makes `x ↦ w(‖x‖)` C² or C⁴ at the centre;
- the first non-zero derivative at `q = 1` has order 4 (C2) or 6 (C4), so the
  zero extension beyond the support is smoother than the centre.

Shape, as chosen by the author: one struct per kernel, `WendlandC2` and
`WendlandC4`, each with a static member template `template <Real T> static
constexpr T operator()(T q) noexcept` (the C++23 static call operator). One
kernel serves `float` and `double`, and the type has no state, which matches
ADR 0003 §4, item 6. The precondition is asserted, as in `Periodic::offset`.

### Numerical requirements

- **Non-negative near `q = 1`.** The expanded polynomial cancels there. For
  C2, `1 − 10q² + 20q³ − 15q⁴ + 4q⁵` returned negative weights at 654 of the
  2,000 largest `float` values below 1, and at 625 of the 2,000 largest
  `double` values. The factored form computes `1 − q` exactly for
  `q ≥ 1/2` (Sterbenz's lemma), so its small values near 1 carry no
  rounding error from that subtraction.
- **Relative accuracy.** Where the exact weight is at least
  `min() / epsilon()` (from `std::numeric_limits<T>`), the relative error is
  a small multiple of `epsilon`. Above that threshold, the factors
  `(1 − q)ᵏ` of the factored forms stay normal. For tiny `q`, terms such as
  `4q` and `c·q` can still be subnormal or round to 0. They enter only sums
  that also contain a term of at least 1, so their absolute error, a small
  multiple of the smallest subnormal, is negligible next to those sums. At
  `q` equal to the smallest subnormal and to `min()`, both kernels return
  exactly 1 in `float` and `double`. Below the threshold, gradual underflow
  may cost relative accuracy, and only the non-negative and finite
  requirement applies.
- **Evaluation.** With `t = 1 − q`:
  - C2: `t2 = t·t`, then `t2 · t2 · (4q + 1)`;
  - C4: `t3 = t·t·t`, then `t3 · t3 · ((c·q + 6)·q + 1)`, with the
    quadratic in Horner form and `c = 35/3` a `constexpr` constant. This
    removes the run-time division by 3, and the polynomial factor is exactly
    1 at `q = 0`. Its coefficients are positive and `q ≥ 0`, so Horner's
    sums do not cancel.

  First-order rounding-error bounds, with `u = epsilon / 2`:

  | Kernel | `q ≥ 1/2` (`t` exact) | `q < 1/2` | Test tolerance `n` |
  | --- | --- | --- | --- |
  | C2 | 5u | 9u | 5 |
  | C4 | 11u | 17u | 9 |

  Each tolerance adds half an `epsilon` for the rounding of the expected
  value. Contracting a multiply and an add into an FMA only removes
  roundings, so the bounds still hold. In a check against exact rational
  values at about 45,000 points, including the 5,000 largest values below 1,
  the largest errors were 3.2 (C2) and 5.4 (C4) `epsilon` in `float`, and 2.2
  and 3.0 in `double`. These are observations, not bounds.
- **Underflow.** At the largest `float` below 1, C4 gives about 8.37e-43, a
  subnormal number. That is allowed: ADR 0003 §4 lets a weight round towards
  0, and the fitter then treats the sample as inactive. In `double` both
  kernels stay normal there.

## K2. The concept

The syntactic requirements, as drafted by the author:

```cpp
template <class K, class T>
concept Kernel = Real<T> && requires(const K &k, T q) {
  { k(q) } noexcept -> std::same_as<T>;
};
```

`T` is a parameter of the concept rather than a member type of the kernel,
so one generic kernel can model `Kernel<K, float>` and `Kernel<K, double>`.

A comment above the concept states the six promises of ADR 0003 §4 that the
compiler cannot check. The concept requires nothing more: no derivative, no
smoothness order, no dimension, no copy requirement and no `constexpr`
check. Add a copy requirement only if Stage 03's fitter stores kernels by
value.

## Tests

`tests/test_kernel.cpp`, with every runtime test run for `float` and
`double`.

### T1. Concept checks

Compile-time checks, for `T` = `float` and `double`:

- Both Stage 01 kernels satisfy `Kernel<K, T>`.
- A `noexcept` lambda on `double` satisfies `Kernel<L, double>` and not
  `Kernel<L, float>`.
- Clean `false` for:
  - `T = int`;
  - a type with no call operator;
  - a call operator that is not `noexcept`, including a lambda without
    `noexcept`;
  - a wrong return type, such as a `float` kernel that returns `double`;
  - a non-`const` call operator.

Fixtures declare their return types. With a deduced `auto` return type,
checking the concept instantiates the body, and an ill-formed body becomes a
hard error instead of a clean `false` (AGENTS.md, P1).

The concept cannot tell `q` from `q²`: a kernel written in `q²` satisfies it.
That rule lives in ADR 0003 §4, item 1, not in a test.

### T2. Kernel values

For each kernel:

- **Ends.** `w(1) == 0` exactly, as ADR 0003 §4 requires, and `w(0) == 1`
  exactly. Both hold for the evaluation above: at `q = 1` the factor `1 − q`
  is exactly 0, and at `q = 0` it is exactly 1 and each polynomial factor is
  exactly 1.
- **Reference values.** A table of `q` values, including `1/2`, whose exact
  weights lie above the relative-accuracy threshold, with the expected
  weights hard-coded as literals. Compute them offline in exact
  rational arithmetic from the expanded polynomial, so the test does not
  repeat the implementation's floating-point expression. Compare with a
  relative tolerance of `n · epsilon`, with `n = 5` for C2 and `n = 9` for
  C4 (see Evaluation), and record the derivation in the test.
- **Non-negative and finite.** On a grid over `[0, 1]`, and on the 2,000
  largest values below 1 (stepping down with `std::nextafter`). This is the
  test the expanded C2 polynomial fails.
- **Positive inside.** `w > 0` on the grid points of `[0, 1)` where the exact
  weight is a normal number.

If the kernels are `constexpr`, `static_assert` checks of `w(0)` and `w(1)`
also confirm compile-time evaluation.

Not tested at runtime, and why:

- **The precondition.** Catch2 has no death tests. The `assert` runs in the
  Debug and sanitizer builds on every in-contract call.
- **Symmetry**, listed in the workflow guide's stage gate. It holds by
  construction: the kernel sees only `q`, so it cannot depend on direction.
- **Smoothness**, also in the stage gate. It is a property of the formula,
  checked symbolically above. The reference values tie the implementation
  to that formula.

## Acceptance

- ADR 0003 is accepted before K1 starts.
- Both kernels satisfy the concept for `float` and `double`, and every
  negative check is a clean `false`.
- T1 and T2 pass with `tools/build-test.sh` (presets `dev-clang`, `dev-gcc`
  and `asan-ubsan`) and in CI.
- `ROADMAP.md` and `docs/INDEX.md` mark Stage 01 complete, and this plan
  becomes reference material, as Stage 00's plans did.

## Open points

- Whether `kernel.hpp` includes all of `field.hpp` only for `Real`, or `Real`
  moves to a smaller header later.
- Whether the concept's parameter becomes `const T q`. In the draft, `q` is a
  non-`const` lvalue, so a kernel taking `T&` satisfies `Kernel<K, T>` and
  could modify the caller's `q`. With `const T q` that kernel is rejected,
  while kernels taking `T` or `const T&` still pass (checked on Clang 18 and
  GCC 13).
