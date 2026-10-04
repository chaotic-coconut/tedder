# Documentation index

Use this table to select relevant documents before reading further.

| Type | Document | Authority/status | Purpose | Useful search terms |
| --- | --- | --- | --- | --- |
| Public README | [README.md](../README.md) | Public project description | Current capabilities, requirements, build commands, and next milestone | C++23, build, periodic, LocalFit |
| Workflow guidance | [AI_DEVELOPMENT.md](AI_DEVELOPMENT.md) | Proposed working method | AI modes, documentation routing, development stages, and delegation boundaries | tutor, reviewer, delegate, stages |
| Accepted ADR | [0001-domain-periodic-contracts.md](adr/0001-domain-periodic-contracts.md) | Accepted decisions; `point_type` decision superseded by ADR 0002 | Contracts and repairs for `field.hpp`, periodic validity, reduction, and documentation | Domain, Periodic, bandwidth, offset, F1, F2, F3, F4 |
| Accepted ADR | [0002-point-type-restriction.md](adr/0002-point-type-restriction.md) | Accepted decisions | Requires `point_type` to equal `Point<scalar_type, dimension>`; supersedes part of ADR 0001 | Domain, point_type, F1 |
| Proposed ADR | [0003-geometry-local-frame-support.md](adr/0003-geometry-local-frame-support.md) | Draft; not authoritative until accepted. Would supersede the `distance_squared` requirement of ADR 0001 F1 | Displacement as the only geometric primitive, scalar bandwidth, internal local frame, support test, neighbour-search contract | Domain, offset, distance_squared, bandwidth, local frame, q², support, search, duplicates |
| Proposed ADR | [0004-fit-results-diagnostics-errors.md](adr/0004-fit-results-diagnostics-errors.md) | Draft; not authoritative until accepted | Derivative order vs degree, `LocalFit` contents, basis, solver and Eigen, diagnostics, `std::expected` errors, non-finite values, float | LocalFit, Hessian, DerivativeOrder, Taylor basis, CPQR, SVD, rank, FitError, std::expected, Clang 19, float |
| Completed plan | [field-plan.md](stages/00-foundation/field-plan.md) | Completed Stage 00 reference material | Production changes and acceptance criteria for the foundation repair | `noexcept`, `is_valid`, minimum image, Jacobian |
| Completed plan | [test-plan.md](stages/00-foundation/test-plan.md) | Completed Stage 00 reference material | Tests, regressions, sanitizer scope, and deferred checks | F1, F2, F3, Catch2, ASan, UBSan |
| Review notes | [2026-10-04-implementation-plan-review.md](reviews/2026-10-04-implementation-plan-review.md) | Not authoritative; input to future ADRs | Review of the 2026-10-04 implementation plan, decisions from discussion, open questions | LocalFit, Hessian, bandwidth limits, distance_squared, support test, std::expected, diagnostics, float |
| Shared instructions | [AGENTS.md](../AGENTS.md) | Project working rules | How AI should work in this repository | learning, delegated, ADR, stage |
| Status map | [ROADMAP.md](../ROADMAP.md) | Current roadmap | Active foundation stage and planned later stages | foundation, kernels, fitting, I/O |
