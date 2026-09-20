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
