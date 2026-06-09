# Stress Testing Singular Gaussian Conditioning Theory

This plan details the implementation of "reverse/stress-test" Julia experiments to validate the theoretical boundaries of the finite-rank shorting theory presented in the paper.

## User Review Required

> [!IMPORTANT]
> The original `bigfloat` validation scripts in `exp2_covariance_calibration` and `exp3_three_node_chain` will remain strictly untouched as the positive control group.

## Proposed Changes

### `experiments/exp4_rough_observation`

#### [NEW] [rough_observation.jl](file:///Users/ostensible_paradox/Documents/math/gaussian-hilbert/experiments/exp4_rough_observation/rough_observation.jl)
- **Objective:** Attack Assumption 4.3 (Cameron-Martin Admissibility). Prove that point evaluation (Dirac delta) fundamentally breaks the theory due to functional-analytic divergence.
- **Implementation:** 
  - Construct a 1D SPDE heat equation prior.
  - Implement a **Narrow Gaussian Approximation** for the observation operator: $g_\varepsilon(x) = \frac{1}{\sqrt{2\pi}\varepsilon} e^{-(x-x_0)^2/(2\varepsilon^2)}$.
  - For a fixed $\varepsilon > 0$, the function is smooth ($g_\varepsilon \in L^2$), satisfying $a_{Ik} \in \operatorname{Range}(Q_0^{1/2})$.
  - Systematically shrink $\varepsilon \downarrow 0$ and compute the exact energy norm.
- **Expected Outcome:** Show that $|Q_0^{-1/2}g_\varepsilon| \to \infty$ as $\varepsilon \downarrow 0$. The $\lambda^{-6}$ remainder estimate breaks down, destroying the asymptotic expansion.
- **Paper Integration:** To be included as **"Example 4: Failure of Cameron-Martin Admissibility"** in the main text.

### `experiments/exp5_constant_diffusion`

#### [NEW] [constant_diffusion.jl](file:///Users/ostensible_paradox/Documents/math/gaussian-hilbert/experiments/exp5_constant_diffusion/constant_diffusion.jl)
- **Objective:** Attack Assumption 4.8 (Dual Scaling). Prove that Fredholm exact cancellation is not a "generic miracle" but relies strictly on the singular clamp scaling mechanism.
- **Implementation:** 
  - Replicate the exact same Lyapunov calibration as `exp2` but fix $D_{\lambda,E} = D_0$ (constant background noise).
  - The drift still scales as $A_{\lambda,E} = -\lambda \Lambda_E$.
- **Expected Outcome:** Show that the Lyapunov equation yields $\Sigma_{EE} \sim \lambda^{-1}$ instead of $\lambda^{-2}$. Consequently, the cross-term $\Sigma_{IE}\Sigma_{EE}^{-1}\Sigma_{EI}$ and the logdet/trace terms no longer share the same asymptotic order. The Fredholm cancellation fundamentally fails.
- **Paper Integration:** Demonstrates that Assumption 4.8 is structurally necessary, not merely a cosmetic assumption. To be promoted to the main text.

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
