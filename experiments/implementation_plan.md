# Stress Testing Singular Gaussian Conditioning Theory

This plan details the implementation of "reverse/stress-test" Julia experiments to validate the theoretical boundaries of the finite-rank shorting theory presented in the paper.

## User Review Required

### `experiments/exp6_galerkin_ill_conditioning`

#### [NEW] [galerkin_ill_conditioning.jl](file:///Users/ostensible_paradox/Documents/math/gaussian-hilbert/experiments/exp6_galerkin_ill_conditioning/galerkin_ill_conditioning.jl)
- **Objective:** Verify Remark 6.1 (previously tracked as Remark 6.3 context): "Why strict Loewner growth is not expected". Show that without strict Loewner monotonicity, the finite-dimensional Galerkin matrix approximations to the Schur complement undergo severe ill-conditioning when adding incremental orthogonal modes.
- **Implementation:**
  - Construct a target observation with a hierarchy of Galerkin subspaces $\mathcal{H}_n$.
  - Compute the finite-dimensional operator residual approximations $R_{\lambda, E}^{(n)}$ as $n$ increases.
  - Track the condition number of the discrete Schur complement matrix $\Sigma_{EE}^{(n)} - R_{\lambda, E}^{(n)}$.
- **Expected Outcome:** Show the condition number of the naive Galerkin Schur blocks exploding, highlighting the necessity of the operator-theoretic inverse-free variational formulation developed in Section 6.2.
- **Paper Integration:** Included as **"Example 6: Galerkin Ill-Conditioning and Loewner Growth"** to justify the abstract variational approach in Theorem 6.1.

## Verification Plan

### Automated Tests
- Run `julia rough_observation.jl` and output error decay curves.
- Run `julia constant_diffusion.jl` and output trace anomaly plots.

### Manual Verification
- Visually inspect the generated `.png` plots to confirm the expected theoretical breakdowns are clearly visible.
