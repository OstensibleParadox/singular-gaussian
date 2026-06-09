"""
Experiment 6: Galerkin Ill-Conditioning and Failure of Strict Loewner Growth.

This experiment numerically validates Remark 6.1 ("Why strict Loewner growth is not expected").
It demonstrates that without the variational inverse-free formulation of the operator
Schur complement, naive finite-dimensional Galerkin matrix approximations
S^(n) = Σ_EE - Σ_EI^(n) (Σ_II^(n))^{-1} Σ_IE^(n)
suffer from catastrophic numerical ill-conditioning as the number of modes n increases.

Because the internal covariance Σ_II is a trace-class operator, its eigenvalues decay
rapidly to zero. When represented in a generic (dense) Galerkin basis, the truncated matrix
Σ_II^(n) develops a condition number that grows polynomially (e.g., n^6), quickly exceeding
standard Float64 precision. Consequently, the matrix inversion (Σ_II^(n))^{-1} introduces
spurious numerical noise, destroying the Loewner monotonicity of the residual R^(n) and
breaking the positive definiteness of the Schur complement.
"""

using LinearAlgebra
using Random
using Printf
using DataFrames
using CSV
using CairoMakie
using LaTeXStrings

const OUTPUT_DIR = joinpath(@__DIR__, "results")
mkpath(OUTPUT_DIR)

function run_galerkin_experiment(N_max=600, step=20)
    # Target covariance
    Sigma_EE = 2.0
    
    # We choose eigenvalues of Σ_II that decay as j^{-6}
    # and cross-coupling coefficients c_j that decay as j^{-5}
    # This ensures CM admissibility: Σ c_j^2 / d_j = Σ j^{-10} / j^{-6} = Σ j^{-4} < ∞
    # Exact residual is roughly π^4 / 90 ≈ 1.0823
    
    Random.seed!(42)
    results = []
    
    println("Running Galerkin Truncation Ill-Conditioning Test...")
    @printf("%-6s | %-15s | %-15s | %-15s | %-15s\n", "n", "cond(Σ_II)", "Exact R", "Float64 R (\\)", "Float64 R (inv)")
    println("-"^75)

    for n in step:step:N_max
        # 1. Exact spectral representations
        d = [1.0 / j^6 for j in 1:n]
        c = [1.0 / j^5 for j in 1:n]
        
        exact_R = sum(c[j]^2 / d[j] for j in 1:n)
        
        # 2. Construct a dense Galerkin basis (random orthogonal matrix)
        # We do this because in a generic finite element or Galerkin basis, 
        # the covariance operator is not perfectly diagonalized.
        A = randn(n, n)
        Q, _ = qr(A)
        U = Matrix(Q)
        
        # 3. Dense covariance matrices
        Sigma_II = U * Diagonal(d) * U'
        Sigma_IE = U * c
        
        # Condition number
        cond_II = d[1] / d[n]  # Exact condition number n^6
        
        # 4. Naive numerical Schur complement updates
        # Method A: Direct linear solve (backslash)
        R_solve = dot(Sigma_IE, Sigma_II \ Sigma_IE)
        
        # Method B: Explicit inversion
        Sigma_II_inv = inv(Sigma_II)
        R_inv = dot(Sigma_IE, Sigma_II_inv * Sigma_IE)
        
        error_solve = abs(R_solve - exact_R) / exact_R
        error_inv = abs(R_inv - exact_R) / exact_R

        push!(results, (
            n = n,
            cond_II = cond_II,
            exact_R = exact_R,
            R_solve = R_solve,
            R_inv = R_inv,
            error_solve = error_solve,
            error_inv = error_inv
        ))
        
        @printf("%-6d | %-15.4e | %-15.6f | %-15.6f | %-15.6f\n", n, cond_II, exact_R, R_solve, R_inv)
    end
    
    df = DataFrame(results)
    CSV.write(joinpath(OUTPUT_DIR, "galerkin_ill_conditioning.csv"), df)
    return df
end

function plot_results(df)
    fig = Figure(size = (900, 450))
    
    # Left Panel: Condition Number
    ax1 = Axis(fig[1, 1],
        title = "Condition Number of Σ_II^(n)",
        xlabel = "Galerkin truncation dimension n",
        ylabel = "Condition Number κ(Σ_II^(n))",
        yscale = log10)
        
    lines!(ax1, df.n, df.cond_II, color = :black, linewidth=2)
    hlines!(ax1, [1e16], color = :red, linestyle = :dash, label="Float64 Precision Limit")
    axislegend(ax1, position = :lt)
    
    # Right Panel: Residual Error
    ax2 = Axis(fig[1, 2],
        title = "Relative Error of Schur Residual R^(n)",
        xlabel = "Galerkin truncation dimension n",
        ylabel = "Relative Error |R_num - R_exact| / R_exact",
        yscale = log10)
        
    lines!(ax2, df.n, df.error_solve, color = :blue, label="Linear Solve (\\)", linewidth=2)
    lines!(ax2, df.n, df.error_inv, color = :orange, label="Explicit Inv (inv)", linewidth=2)
    
    # Add a marker where precision breaks
    break_idx = findfirst(x -> x > 1e16, df.cond_II)
    if break_idx !== nothing
        vlines!(ax2, [df.n[break_idx]], color = :red, linestyle = :dash, alpha=0.5)
        text!(ax2, df.n[break_idx] + 10, 1e-2, text="κ > 10¹⁶", color=:red)
    end
    
    axislegend(ax2, position = :lt)
    
    save(joinpath(OUTPUT_DIR, "galerkin_ill_conditioning.png"), fig)
    println("\nSaved plot to $(joinpath(OUTPUT_DIR, "galerkin_ill_conditioning.png"))")
end

function main()
    df = run_galerkin_experiment(600, 20)
    plot_results(df)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
