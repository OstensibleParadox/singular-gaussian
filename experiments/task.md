# Stress Tests Implementation Tasks

- [x] Create 3 independent experiment directories (`exp4_rough_observation`, `exp5_constant_diffusion`, `exp6_galerkin_ill_conditioning`).
- [x] Implement Exp 4 (Failure of Cameron-Martin Admissibility).
  - [x] Write `rough_observation.jl` using the narrow Gaussian approximation $g_\varepsilon(x)$.
  - [x] Compute and plot the divergence of the Cameron-Martin energy norm $|Q_0^{-1/2}g_\varepsilon| \to \infty$ as $\varepsilon \downarrow 0$.
  - [x] Show the breakdown of the $\lambda^{-6}$ remainder estimate.
- [x] Implement Exp 5 (Failure of Dual Scaling / Constant Diffusion).
  - [x] Write `constant_diffusion.jl` with $D_\lambda = D_0$.
  - [x] Verify $\Sigma_{EE} \sim \lambda^{-1}$ scaling under the Lyapunov equation.
  - [x] Generate trace anomaly plots showing the asymptotic decoupling and failure of Fredholm exact cancellation.
- [x] Implement Exp 6 (Galerkin Ill-Conditioning).
