# RustSolid modeling operator contracts (v0.12)

No external CAD kernel is used as a mathematical oracle. A passing `Solid::validate()` is necessary but does not prove that an operator produces the intended geometry.

The executable inventory is `src/operator_contract.rs::ALL_MODELING_OPERATORS`; the pinned witnesses are in `tests/fixtures/operator_contracts.v1.json` and their implementations in `tests/operator_contracts.rs`. Each registered operation needs at least one independent analytic test, one metamorphic test, and two invalid-boundary cases.

The current gate covers three model constructors (`block`, `extrude_xz`, `analytic_cylinder`) and all six `EulerCommand` variants. Primitive cases refer to permanent reviewed analytical fixtures in `tests/fixtures/regressions.v1.json`. Euler cases execute the actual Command/Memento pathway on a cuboid and verify the closed-form volume, area, centroid, bounds, Euler characteristic, expected topological changes, Undo and Redo. Metamorphic scenarios verify translations and the square/cube scaling laws. Invalid cases verify that failures do not mutate geometry or history.

**Enforcement:** adding a new `EulerCommand` variant without extending its exhaustive `operator_kind()` mapping causes a Rust compilation error. Updating the Rust operation inventory without the corresponding versioned test scenarios fails the operator contract CI gate. The test also deliberately removes required categories and substitutes forged witnesses to establish that the gate rejects incomplete coverage. New public primitive constructors must be explicitly entered in the inventory; Rust does not automatically reflect arbitrary public function names.

## Adding an operation

1. Add its typed `ModelingOperator` ID and register it in `ALL_MODELING_OPERATORS`; for commands, extend `EulerCommand::operator_kind()`.
2. Add at least four executable scenarios to `operator_contracts.v1.json`. Verify independent expected values mathematically, not against another CAD implementation.
3. Write a real scenario driver in `tests/operator_contracts.rs`; adding labels without executable behavior is invalid.
4. Add boundary cases, include tolerances and scales and check that invalid operations are atomic. Promote every newly discovered defect to `regressions.v1.json` with exact inputs.
5. Run `cargo test --test operator_contracts` and `cargo test --all-targets`; CI enforces both.

The gate covers its listed operators and test families only. It is not a mathematical proof of industrial robustness or a substitute for deep testing of arbitrary Boolean, NURBS and trimmed-surface operations.
