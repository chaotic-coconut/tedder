# tedder:

Revised 2026-09-19. This replaces the
previous testing plan. Read [FIELD_PLAN.md](FIELD_PLAN.md) first: it defines the
production contracts. These documents are instructions, not implemented changes.

## Scope and order

Keep this repair small. Add the three tests in A, implement F1–F4 together with
their regressions in B, then add the sanitizer job in C. Preserve the existing
tests. Do not build a general testing framework or start reconstruction work.

The reviewed baseline was
[c11b602](https://github.com/chaotic-coconut/tedder/tree/c11b60236c0d4634fefd5f72fe75b235e1f81fa6).
Inspect the actual checkout before editing and retain subsequent fixes. Use the
existing Catch2 setup; one additional `test_types.cpp` is enough for A. Register
any new sources in `tests/CMakeLists.txt`.

| Tests | Production dependency |
| --- | --- |
| A: matrix, result layout, aligned subview | None; should pass before repairs |
| B1: tightened `Domain` | F1 |
| B2: invalid periodic configuration | F2 |
| B3: periodic reduction regressions | F3 |
| Documentation review | F4 |

Do not make tests preserve the known bugs. If implementation reveals a contract
conflict, resolve it in both plans before encoding a different expectation.
For a tests-only assignment, complete independent work and report F1–F3 tests as
pending their corresponding repairs; do not silently change production contracts.

## A. Three tests of existing code

1. **Rectangular matrix:** use `Matrix<double,2,3>`. Set six distinct entries
   through `operator[]` and compare `data()` with an explicit row-major sequence.
   Read through a const reference as well. Do not compute expected indices using
   the same expression as the implementation.
2. **Rectangular Jacobian:** use `LocalFit<double,2,3>` (`D=2`, `C=3`). Assert that
   its Jacobian has three rows and two columns, and write/read `[2,1]`. This catches
   a swapped `D,C` declaration that square examples miss. Do not invent fitting
   semantics for `neighbors` or `bandwidth`.
3. **Aligned subview:** construct a `SampleView` from four distinct point/value
   pairs in live owning arrays. Take `subview(1,2)` and verify both selected pairs,
   its size, and the starting addresses of both spans. This checks tedder's
   selection and alignment; an owner-mutation demonstration is unnecessary.

No exhaustive cross-product of types, dimensions, and subview shapes is required.

## B1. F1: Domain constraints

Keep the existing positive checks for both built-in models. Add a small valid
custom domain whose `point_type` is a distinct point struct, not `tedder::Point`.
It must satisfy `Domain`; point representation remains deliberately unrestricted.

Using otherwise-valid fixtures, require clean `false` results for:

- an integral scalar and a `std::string` scalar;
- a mutable runtime `dimension`, zero dimension, and negative signed dimension;
- each of `offset`, `distance_squared`, and `max_bandwidth` lacking `noexcept`;
- an unrelated type such as `int`.

Preserve the existing missing-method, wrong-return-type, and non-const fixtures.
Add `noexcept` to their otherwise-valid methods so these checks still isolate
the intended defect. A fixture rejected for two unrelated reasons can hide a
broken constraint. Give malformed-scalar fixtures their own point type rather
than instantiating the separately constrained `Point` alias with an invalid type.

Use ordinary `static_assert(!Domain<Bad>)` checks: these must compile without hard
template errors. This is not the deferred compile-fail harness. Confirm that the
valid custom domain works with the free helpers and preserves their `noexcept`
contract. Do not execute a throwing-domain termination test.

## B2. F2: periodic validity and bandwidth

Check valid default, finite mixed, and all-open configurations. Finite negative
periods still mean open axes. For NaN, positive infinity, and negative infinity
in `length`, require all three outcomes:

- `is_valid()` is false;
- `max_bandwidth()` is zero;
- `admits_bandwidth()` rejects the tested bandwidths.

Use a small table, including an invalid entry after a valid positive period, so
validation cannot accidentally stop at the first usable axis. Exercise mutation
from valid to invalid and back; also change a positive period and check that the
bandwidth limit updates. Nothing may cache the old public `length` contents.

For a valid period of 2, check `nextafter(1,0)`, 1, and
`nextafter(1,+infinity)`: only the first is admitted. Retain rejection of zero,
negative values, NaN, and both infinities. On an invalid domain, include a normal
positive bandwidth among the rejected inputs.

Do not call `offset` or `distance_squared` on an invalid configuration in the
ordinary suite. That violates a documented precondition, not an exception
contract. The assertion subprocess test remains deferred.

## B3. F3: periodic reduction

Add four named regression cases. The inputs below target the binary `float` and
`double` formats used by the current CI toolchains.

| Case | Type | Period `L` | `p` | `q` | Expected offset |
| --- | --- | --- | --- | --- | --- |
| Large double coordinate | `double` | 3 | `1e16` | 0 | 1 |
| Large float coordinate | `float` | 3 | `1e8f` | 0 | 1 |
| Rejected wide fast-path threshold | `float` | `nextafter(1.0f,2.0f)` | `4000.25f` | 0 | `0.249523162841796875f` |
| Overflowing raw subtraction | `double` | `1e308` | `1e308` | `-1e308` | 0 |

The third expected value is independently derived as
`1/4 - 4000/2^23`, exactly representable in binary32. The fourth uses opposite
multiples of the same stored period. Check the offset directly; it does not
establish general correctness of squared distances at enormous magnitudes.

Check reversed arguments for these cases. Keep the existing negative-coordinate,
outside-box, mixed-axis, offset-sign, and canonical half-period tests. For `L=2`,
the existing convention remains `offset({0},{1}) == {1}`, with reverse `{-1}`.

Add a small set of representable examples around the fast-path boundary
`abs(p-q) == L` and the half-period boundary. Cover both signs and both branches,
using `nextafter` where useful. Check finite bounded offsets and antisymmetry.
Test period-shift invariance away from half-period ties; at a tie, equivalent
coordinate representations may select opposite equally valid minimum images.
Do not require global representation-invariant offset signs there.

Use exact equality for the four regressions and deliberately exact examples.
Elsewhere use an absolute tolerance tied to the scalar type and period, such as
`8 * epsilon(T) * L` for ordinary single-axis examples. Include an absolute
tolerance near zero. Do not use one double-only tolerance for float tests or
loosen a tolerance to accept the old wide-threshold error.

Expected values must come from known geometry or independently derived
remainders. Do not copy the production reduction into a test oracle. A large
random sweep or benchmark is not an acceptance gate for this repair.

## C. Build and sanitizers

Retain the existing GCC 13 / Clang 18 Debug and Release jobs. Add one separate
Linux Debug job running the complete suite with ASan and UBSan:

- instrument compilation and linking with `-fsanitize=address,undefined`;
- retain stack frames with `-fno-omit-frame-pointer`;
- make UBSan findings fail the job, for example with
  `-fno-sanitize-recover=undefined`;
- run CTest with `--output-on-failure --no-tests=error`.

Keep sanitizer flags local to the diagnostic build/test target; do not export
them through `tedder::tedder`. Use the existing CMake workflow, explicitly enabling
`TEDDER_BUILD_TESTS`. No consumer project, additional compiler matrix, or CI
redesign is required. Catch2 assertions must remain active in Release.

Run the configurations available locally and let CI cover the remainder. Report
unavailable runs as unverified. Sanitizers cover executed memory and undefined
behaviour errors; they do not verify wrapping accuracy or prove lifetime safety.

## Deferred checks and limits

| Check | Why it is not required now |
| --- | --- |
| Exhaustive matrix/view coverage and `long double` matrix | Beyond the agreed small regression scope |
| Assertion subprocesses | Possible, but need a separate process; `assert` does not throw and Release need not reject misuse |
| Hard compile-fail tests such as zero matrix extents | Possible, but need a compiler harness that forces instantiation; ordinary negative concept checks remain required |
| Consumer CMake project and two-translation-unit linkage | Useful integration checks, deferred from this repair |
| Dangling views and out-of-bounds access in ordinary tests | Invalid use has no supported output; a passing sanitizer run cannot prove every caller's storage lifetime |
| Periodic constant evaluation on every toolchain | Compiler and standard-library support varies; preserve existing supported checks and runtime coverage, without adding a capability-check framework now |
| Exhaustive numerical guarantees | Finite examples cannot prove all magnitudes and sample layouts; documented range limits still apply |

Do not replace deferred tests with empty passing placeholders. `std::sqrt` is not
portably constexpr in C++23; runtime distance coverage is sufficient here. See
F4 for the standard/library distinction.

Fitting tests wait for fitting code. When that work starts, follow `PROJECT.md`:
fixed bandwidth first, polynomial reproduction within tolerance, analytic
derivatives with physical scaling, and explicit failure for insufficient or
degenerate neighbourhoods. Do not add solver, kernel, kd-tree, or I/O work here.

## Completion report

Report files changed, F1–F4 status, the regressions added, and the exact build/CI
configurations actually run. Identify remaining failures and unavailable checks.
Review F4's comments against the implementation; passing tests alone do not prove
the documented contracts.
