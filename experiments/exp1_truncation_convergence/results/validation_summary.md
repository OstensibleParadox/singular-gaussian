# Experiment C.1 Validation and Truncation Convergence Summary

## Claim Tested

This Julia experiment evaluates truncation convergence for the trace-class OU witness.

## Truncation Diagnostics

| Quantity | Value |
|---|---:|
| Reference Schur scalar at lambda=4 | 3.066990926746e-02 |
| m2 extrapolation at n=1000 (Float64) | 4.999999999999e-01 |
| Absolute error to epsilon0/2 at n=1000 | 6.867e-14 |
| Baseline Lyapunov max residual at N_ref | 5.551e-17 |
| Intervened Lyapunov max residual at N_ref | 1.388e-17 |

| Diagnostic slope over n >= 128 | Value |
|---|---:|
| Noise trace tail | -0.499 |
| Baseline covariance trace error | -2.683 |
| Intervened covariance trace error | -2.683 |
| Schur scalar error | -5.504 |


## Output File Checksums

| File | MD5 | SHA-256 |
|---|---|---|
| truncation_convergence.csv | 2a7ef7db3ca7ebb1dce77c2663fe6725 | c56d052c3342bb488bde93a5c7796d8b194745b2d694e0208dbcc515e5ae8b14 |
| lambda_schur_scaling.csv | 2a8beaad700a5be9fc777c6ed02c177b | e0309eabd58d55a5633da04563f3844bdf054cdd9b1da98e20f46548541337cd |
| truncation_convergence_validation.png | d894ac7eb4ff79ac0631c29a3b2fff63 | a2b9552c2cff9deede913a4fe53d40118b347712b0b97ddc233a98544144ead8 |
