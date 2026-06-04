# Experiment C.2 validation summary

The experiment evaluates the target covariance calibration ratio `lambda * d_lambda / eps0`
and the defect `|lambda * d_lambda / eps0 - 1|`.

Core evaluation precision: `256` bits (approximately `77` decimal digits).

| family | ratio(lambda=1) | ratio(lambda=4096) | defect(lambda=4096) | result |
|---|---:|---:|---:|---|
| Canonical | 1 | 1 | 0 | exactly scale invariant |
| Slow decay | 1 | 4096 | 4095 | ratio diverges |
| Fast decay | 1 | 0.00024414062 | 0.99975586 | ratio collapses to zero |
| Asymptotic perturbation | 2 | 1.0000001 | 5.9604645e-08 | asymptotically correct, not scale invariant |

The asymptotic perturbation family reaches the calibrated limit but fails the
exact finite-lambda invariance condition. The experiment checks the consequence
of the covariance calibration condition.

High-precision probe: at `lambda = 2^100`, the perturbed
family remains non-constant with exact-calibration defect = 6.223015277861141707144064053780124240590252168721167133101116614789698834035383e-61.

## Output File Checksums

| File | MD5 | SHA-256 |
|---|---|---|
| ablation_data.csv | b5057f495d190966e39091cf5076603d | ac68ecc7202d3546fbbdc1282f69d7091163f13a4287dba7bca415db6032a869 |
| precision_probe.csv | a4017c6eb084bf24f1387e169a841731 | e8bc3395e9134ee5f4ab34c316b0bdabd4e9c0b5bdc10bda3149f9f419427fac |
| ablation_combined.png | 192ea099348c18765630af69cb91eba9 | 0095c73b6cecf5a51be4a04afc6d52a0ce0cb2a66de922b5c8e75603895aa900 |
