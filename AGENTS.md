# tedder agent instructions

- tedder is a modern-C++ learning project.
- The user normally designs and implements learning-heavy core library code.
- For core code, AI defaults to tutor or reviewer mode.
- AI may implement mechanical work only when explicitly delegated.
- Accepted ADRs override broader design documents.
- Active stage plans define current implementation work.
- Do not turn undecided policies into public APIs.
- Before non-trivial work, read [docs/INDEX.md](docs/INDEX.md) and then only
  the relevant routed documents; do not automatically load every project document.
- Work on one roadmap stage at a time and preserve unrelated user changes.
- Keep changes narrow and do not silently rewrite accepted decisions.
- Run the smallest relevant checks while working and report commands actually run.
- Do not claim numerical identities or guarantees without support from the
  applicable contract and representation.
  
## Review guidelines

- Check changes against accepted ADRs in `docs/adr/` and the active stage plan.
  A contradiction with an accepted ADR is P1. Proposed ADRs are context, not rules.
- Treat these as P1:
  - a claim of an exact floating-point result, identity or guarantee that the
    contract and representation do not support;
  - overflow, underflow, NaN or tie behaviour that silently changes a result;
  - a tolerance or rank threshold written as a fixed literal instead of a
    multiple of `std::numeric_limits<T>::epsilon()`;
  - a test whose expected value is computed with the same formula as the code
    under test;
  - a public API change (concept requirements, public types, signatures,
    `noexcept`, `constexpr`) without an ADR or stage-plan entry;
  - a negative concept check that would be a hard compile error rather than a
    clean `false`.
- Every behaviour change needs tests; code templated on the scalar type needs
  tests for both `float` and `double`.
- Ignore pure formatting; `.clang-format` decides it.
- This is a learning project: explain each finding and why it matters, and do
  not propose wholesale rewrites of core code.
