# tedder: Stage 02 neighbour-search plan

Status: active. Decisions A1 to A5 are made; S1 can start.

This plan implements the parts of
[ADR 0003](../../adr/0003-geometry-local-frame-support.md) that turn a query
point and a bandwidth into a candidate set:

- §3: the isotropic local frame, forward direction only (`r ↦ ξ`);
- §4: the support variable `q²` and the candidate test `q² < 1`;
- §5: the search contract, a brute-force backend, and the step that turns
  backend output into candidates (sort, remove duplicates, support test,
  count).

The ADR is authoritative; this plan does not repeat its reasons.

## Scope

In scope:

- the frame, `q²` and the candidate test;
- a brute-force backend that meets the §5 contract;
- candidate selection: sort, remove duplicates, test, count;
- their tests in `tests/test_search.cpp`, registered in
  `tests/CMakeLists.txt`, for `float` and `double`.

Out of scope:

- kernel weights, the active-sample count, and fitting (Stage 03);
- transforming derivatives back to physical coordinates (Stage 03, §3);
- a public search concept. It waits for the second backend, nanoflann, in
  Stage 05, as the workflow guide asks;
- frames other than `B = I/h`, per-axis or variable bandwidths;
- validating `h` and `x₀` with reported errors. That is the fit's job in
  Stage 03 (ADR 0004 §7); here they are preconditions.

`field.hpp` and `kernel.hpp` do not change.

## Decisions

| ID | Question | Decision |
| --- | --- | --- |
| A1 | Brute force: return every index, or prefilter with the canonical `q² < 1` test? This was ADR 0003's open point. | Every index. In the author's words, brute force is only for regular grids and maybe a fallback. It is complete without any argument, it does not depend on the frame, and candidate selection stays the one place that decides. |
| A2 | Search output: a caller-owned buffer, or a returned vector? | A caller-owned `std::vector<std::size_t>&`, reused across queries. The backend clears and refills it, which keeps its capacity, so a query that returns no more indices than the buffer already holds does not allocate. |
| A3 | Public surface in Stage 02. | None: everything in `namespace detail`. The fit in Stage 03 is the first public entry point, and the search concept comes with Stage 05. |
| A4 | Spelling and names. | American spelling in code: identifiers, header and file names, as in the C++ standard library (`initializer_list`), Eigen (`normalize`) and nanoflann (`findNeighbors`). British prose stays as it is, and `docs/INDEX.md` lists both spellings as search terms. The rule is also in `AGENTS.md`. The proposed names (`IsotropicFrame`, `BruteForceSearch`, `detail/frame.hpp`, `detail/search.hpp`, `tests/test_search.cpp`) avoid the word either way. |
| A5 | Non-finite sample coordinates. | A precondition, as in ADR 0001, with a debug-only assertion in candidate selection. Whether it becomes a reported failure is decided in Stage 03, where errors are reported. |

## Order of work

| Step | Work | Tests |
| --- | --- | --- |
| S1 | Frame and `q²` | T1 |
| S2 | Candidate test | T2 |
| S3 | Candidate selection | T3 |
| S4 | Brute-force backend | T4 |
| S5 | End-to-end configurations | T5, T6 |

The author designs and writes S1 to S4. AI may write the tests when
delegated, and reviews units and inclusivity (whether the boundary is in or
out), as the workflow guide's stage table says.

## S1. Frame and `q²`

- **Frame.** A concrete `detail` type for `B = I/h`. Its only Stage 02
  operation, `to_local`, maps a physical offset `r` to `ξ`, with
  `ξ_d = r_d / h`: one division per component, as ADR 0003 §4 fixes. It
  never multiplies by `1/h`.
- **Names.** `to_local` names the role, not the arithmetic, so it stays
  right for a general `B` with rotation. Stage 03 adds the reverse transforms
  for Jacobians and Hessians as `to_physical`. The argument is an offset,
  not a point: `offset` has already translated, and the frame is linear
  (ADR 0003 §3).
- **Type name.** `IsotropicFrame`, decided by the author. Offsets form a
  vector space, whose frame is a basis with no origin, so the type needs no
  `x₀`. `Isotropic` separates it from a possible general-`B` type, which
  could be `MatrixFrame`, with `LocalFrame` for a later `detail` concept
  (ADR 0003 §3). Not chosen: `IsotropicScaling`, accurate but leaving no
  family name for the concept; `IsotropicBandwidth`, which names the
  parameter, while in kernel density estimation a bandwidth matrix is a
  variance matrix (`H = h²I`), not `B`; `IsotropicNormalization`, because
  ADR 0003 §4 already uses "normalising" for the kernel's constant.
- **`q²`.** `q² = Σ ξ_d²`, summed in axis order `0 … D−1`.
  - It is not `distance_squared(g, x, x₀) / h²`. That is the rejected
    `r · r < h²` form, which fails for extreme `h`.
- **Exactness, for tests only.** `r_d / h` is exact when `h` is a power of
  two and the result stays normal.
- **Precondition.** `admits_bandwidth(g, h)`, asserted where the frame is
  built.

Proposed shape; the names `IsotropicFrame` and `to_local` are decided:

```cpp
namespace tedder::detail
{
  template <Real T, std::size_t D>
  struct IsotropicFrame
  {
    T h;
    constexpr Point<T, D> to_local(const Point<T, D> &offset) const noexcept; // ξ = offset / h
  };

  template <Real T, std::size_t D>
  constexpr T squared_norm(const Point<T, D> &xi) noexcept; // q², axis order
}
```

## S2. Candidate test

A sample is a candidate exactly when `q² < 1`, written so that a NaN `q²`
fails it (ADR 0003 §4). The boundary is excluded: `q² = 1` is not a
candidate. The same function is used everywhere a candidate is decided.

## S3. Candidate selection

Input: the domain, the frame, `x₀`, the sample points, and the backend's
index buffer. In place, it:

1. sorts the indices and removes duplicates (`std::sort`, `std::unique`,
   `erase`);
2. keeps the indices whose `q²` passes S2, with `r = g.offset(x_i, x₀)`.

The buffer then holds the candidates in ascending sample index, and its size
is the candidate count of ADR 0004 §6.

- **Preconditions:** every index is below the number of samples, and every
  sample it reads has finite coordinates (both asserted, A5); `x₀` is
  finite.
- **Kept for Stage 03:** `ξ` and `q²` are not stored. Stage 03 decides
  whether to recompute them or to make candidate selection return them
  alongside the indices. Recomputing uses the same functions, but gives the
  same `q²` only if both calls are compiled alike (see Open points, FMA
  contraction).

## S4. Brute-force backend

The reference backend for Stage 05's equivalence tests. It meets the §5
contract: complete, extras allowed, duplicates allowed, no decisions.

With A1, a query replaces the buffer's contents with
`0, 1, …, N−1`. It holds the domain and a `PointView` of the samples. The
storage the view refers to must outlive the backend, as for `SampleView`.

Proposed shape:

```cpp
namespace tedder::detail
{
  template <Domain G>
  class BruteForceSearch
  {
  public:
    BruteForceSearch(const G &g, PointView<typename G::scalar_type, G::dimension> points);
    void query(const typename G::point_type &x0, typename G::scalar_type h,
               std::vector<std::size_t> &indices) const;
  };
}
```

Indices are `std::size_t` positions in the `PointView`. A Stage 05 adapter
converts nanoflann's index type.

## Tests

`tests/test_search.cpp`. Every test runs for
`float` and `double`. No expected value is computed with the formula under
test: expectations are hand-built, or come from integer arithmetic.

### T1. Frame and `q²`

- With `h` a power of two, `ξ` equals hand-computed dyadic values exactly,
  including negative components.
- `to_local` divides. With `h = 3`, `5 / h` equals `5/3` rounded once, from
  exact arithmetic. Multiplying by a rounded `1/h` differs in the last bit,
  in `float` and in `double`.
- `q²` of hand-chosen `ξ` with exactly representable squares and sums.
- `q²` sums in axis order. With `ξ = (1, a, a, a, a)` and `a²` at most half
  an ulp of 1 (`a = 2⁻¹²` in `float`, `2⁻²⁷` in `double`), the sum in axis
  order is exactly 1; summing the `a²` first gives more.

### T2. Candidate test

| Case | Expected |
| --- | --- |
| `r = (h, 0, 0)`, `h` a power of two: `q² = 1` exactly | not a candidate |
| `r = (nextafter(h, 0), 0, 0)` with the same power-of-two `h`: `ξ₀` is the largest value below 1 and `q² < 1` | candidate |
| `q² = 0`, `r = 0` | candidate |
| `q²` NaN, passed to the test directly | not a candidate |
| `q² = +∞` | not a candidate |
| Far outside: `h = min()`, `r = 1`. `q²` overflows to `+∞` | not a candidate |
| Far inside: `h = max()`, `r = min()`. `ξ` underflows to 0 | candidate |

The predecessor row needs `h` to be a power of two. For other `h` the
quotient can land two steps below 1: with `float` `h = nextafter(1, 2)` it is
`0x1.fffffcp-1`.

The NaN case calls the test with a NaN `q²` directly. It does not pass
non-finite coordinates, which would break a precondition. The last two rows
are the range examples of ADR 0003 §4; both behave as listed in `float` and
in `double`.

### T3. Candidate selection

Hand-built index buffers with a hand-listed expected result:

- unsorted indices with duplicates give sorted, unique candidates;
- indices outside the support (extras) are removed;
- an empty buffer gives an empty result;
- two samples at the same coordinates are both kept, as distinct indices;
- the sample at `x₀` itself is kept.

### T4. Brute-force backend

For `N = 0`, `1` and a few more, a query returns every index exactly once,
whatever `x₀` and `h`. The buffer's previous contents are replaced.

A second query of the same size leaves the buffer's `data()` and
`capacity()` unchanged (A2). This catches an implementation that builds a new
vector and assigns it. It cannot see one that frees and reallocates at the
same address; a counting `operator new` would, but it replaces the global
allocator for the whole test executable, so it is not planned.

### T5. Hand-built configurations

The backend followed by candidate selection, against hand-listed candidate
indices:

- `Euclidean` in 1, 2 and 3 dimensions;
- **periodic seam:** `Periodic` with `L = 1`, `x₀ = 1/32`. A sample at
  `31/32` is at offset `−1/16`, a candidate for `h = 1/8`. The mirror case is
  `x₀` near 1;
- **mixed domain:** periodic in `x`, open in `y`, with a sample that is near
  only across the seam in `x`;
- **half-period tie:** `L = 1`, `h` the largest value below `1/2`, and a
  sample exactly `1/2` away. It is not a candidate (ADR 0003 §4);
- **no candidates:** samples present, none in the support;
- **no samples;**
- **extreme bandwidth:** `Euclidean`, `x₀ = 0`, `h = min()`, samples at `0`,
  `h/2` and `2h`. The first two are candidates and the third is not. The
  rejected forms `r · r < h²` and `distance_squared / h²` reject the first two,
  because `h²` underflows to 0, so this case checks that selection uses the
  canonical `ξ` (ADR 0003 §4).

Except in the extreme-bandwidth case, all coordinates are multiples of `1/32`
in `[0, 1]`, and `L = 1`. The difference of two such values is a multiple of
`1/32` no larger than 1 in magnitude, so it is representable and the
subtraction is exact; so are the periodic reductions. Being dyadic alone would
not be enough: in `float`, `2^(−25) − 1` rounds to `−1`.

### T6. Dyadic grid against an integer oracle

- **Grid:** `h = 2^(−3)`, `x₀` a multiple of `h` in `[0, 1]`, and samples
  at `x₀ + m · 2^(−j) · h` for every integer vector `m` with
  `|m_d| ≤ 2^(j+1)`. With `j = 3` in 3 dimensions that is `33³ = 35,937`
  samples per `x₀`, checked exhaustively. Every coordinate, offset, `ξ`,
  square and partial sum is exact in `float` and in `double`.
- **Oracle:** the expected decision is integer arithmetic,
  `Σ m_d² < 4^j`. It is independent of the code under test, and it includes
  the six boundary points, `m = (±8, 0, 0)` and its permutations.
- **Domains:** run it for `Euclidean`, and for `Periodic` with `L = 1`,
  `x₀` near 0 and the sample coordinates wrapped into `[0, 1)`, so that the
  offsets cross the seam.

Checked exhaustively for three values of `x₀`, for `Euclidean` and for
`Periodic` with wrapped coordinates, in both types: every offset was exact
and every decision matched the integer oracle.

Points off such a grid have no exact boundary: ADR 0003 §4 records
misclassification within about 3 ulps of `q = 1` for every formula. So no test
compares against exact arithmetic there.

## Acceptance

- The decisions A1 to A5 are recorded here (done before S1).
- T1 to T6 are in `tests/test_search.cpp`, registered in
  `tests/CMakeLists.txt`, and pass with `tools/build-test.sh` (presets
  `dev-clang`, `dev-gcc` and `asan-ubsan`) and in CI.
- No public API is added.
- On completion:
  - `ROADMAP.md` and `docs/INDEX.md` mark Stage 02 complete;
  - `README.md` updates its status and its "Next" line;
  - ADR 0003's Consequences record what Stage 02 implemented, and keep only
    Stage 03's part of §3 to §5 as still to do;
  - this plan becomes reference material.

## Open points

- Whether candidate selection should return `ξ` and `q²` with the indices.
  This is decided in Stage 03, by what the fit needs.
- FMA contraction in `squared_norm`. A compiler may fuse
  `sum += xi[d] * xi[d]` into a fused multiply-add, depending on flags and
  target. clang-18 for ARM64 does by default, and x86-64 builds do with
  `-march=haswell`. Default x86-64 builds, as in CI, do not, and neither
  does constant evaluation with gcc-13 or clang-18. So `q²` can differ
  between platforms, and between compile time and run time.
  - The difference grows with `D`. Without fusion, each of the `D` products
    rounds once more, and the roundings accumulate. Both sums lie within
    `γ_D · q²` of the exact `q²`, with `γ_D = D·u / (1 − D·u)` and
    `u = epsilon / 2`. So two builds can decide differently only when the
    exact `q²` is within about `D·u` of 1. Example from the review of PR
    #12: `float`, `D = 22`, every `ξ_d = 0x1.b4a29p-3`. The unfused sum is
    `0x1.ffffecp-1` and the fused sum exactly 1, although the exact `q²` is
    about `4.4 u` below 1.
  - ADR 0003 §4 accepts that no formula decides exactly near the boundary.
    Its "about 3 ulps" is an observation, not a bound. What it requires is
    that every caller uses one formula. Contraction lets that one formula
    round differently in different builds, so the question is consistency,
    not accuracy.
  - The planned tests use exact values, where fusing changes nothing.
  - The rest of the path is unaffected. Division is correctly rounded, and
    `Periodic::offset` multiplies `L` by −1, 0 or 1, which is exact.
  - If Stage 03 recomputes `q²` and needs the same bits as selection, an
    explicit `std::fma` in `squared_norm` fixes its rounding on every
    platform with IEEE arithmetic, that is, without `-ffast-math`.
    tedder is header-only, so it cannot impose `-ffp-contract=off` on its
    users' builds. The costs: a call to the math library per term in
    default x86-64 builds, and no compile-time evaluation with clang-18.
  - If candidate selection returns `q²`, nothing is recomputed and the
    question disappears. Decided together with the point above.
- Whether a per-sample stamp array should replace sort and unique for very
  large clouds (review D9). Not before measurements show the sort matters.
