# Experiment C.3 validation summary

The experiment evaluates the three-node Gaussian SDE X_1 -> X_2 -> X_3.
It validates the target Schur residual scaling lambda^2 * S_{2, lambda} -> 1/2
and the renormalized clamped action J_lambda under three cases.

Core evaluation precision: `256` bits.

| Case | description | lambda^2 * S_2(lambda=4096) | J(lambda=4096) | J(lambda=2) |
|---|---|---:|---:|---:|
| A | Decoupled-centered | 0.50000000 | 0 | 0.00000000 |
| B | Coupled-centered | 0.50000000 | 1.77548979517e-15 | 0.01298774 |
| C | Coupled-off-centered | 0.50000000 | 1 | 1.01298774 |

Results show:
- Panel A: lambda^2 * S_{2, lambda} converges to 0.5 in all cases.
- Panel B: J_lambda is exactly 0 under Case A, decays as O(lambda^-4) to 0 under Case B,
  and converges to 1.0 under Case C.

## Output File Checksums

| File | MD5 | SHA-256 |
|---|---|---|
| three_node_data.csv | 8b7160f3e1ced09d89f44817bd8af5be | a7cfd4a05ce97070c038c3a78ae12c2d7266c9268f63e67af4182e782f561963 |
| three_node_gaussian_chain.png | 5eafea4c8a98d5b30c0511fe7e9b2e94 | c43061d633a4bffe4851b3a8f63f09008eaba4bae5d3be7f3199aeb36a0bd6aa |
