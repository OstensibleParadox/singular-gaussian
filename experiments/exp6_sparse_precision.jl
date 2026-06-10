using LinearAlgebra
using SparseArrays
using Printf
using Pkg

# Ensure Plots is installed
if !("Plots" in keys(Pkg.project().dependencies))
    Pkg.add("Plots")
end
using Plots

function build_laplacian(N, kappa, dx)
    diag_elements = fill(2.0 / dx^2 + kappa^2, N)
    off_diag_elements = fill(-1.0 / dx^2, N - 1)
    L = spdiagm(-1 => off_diag_elements, 0 => diag_elements, 1 => off_diag_elements)
    return L
end

function get_fourier_vector(N, dx, exponent)
    v = zeros(N)
    for j in 1:N
        if j % 2 != 0
            v[j] = 1.0 / (j^exponent)
        end
    end
    v_spatial = zeros(N)
    for i in 1:N
        x_i = i * dx
        val = 0.0
        for j in 1:N
            if v[j] != 0
                val += v[j] * sin(j * pi * x_i)
            end
        end
        v_spatial[i] = val
    end
    return v_spatial ./ norm(v_spatial)
end

function run_experiment_and_plot()
    N = 10000
    kappa = 1.0
    Delta_E = 1.0
    dx = 1.0 / (N + 1)
    
    L = build_laplacian(N, kappa, dx)
    
    # We use the Rough Observation (1/j^3) which violates Cameron-Martin for GFF
    # but satisfies finite energy for OU.
    v = get_fourier_vector(N, dx, 3)
    Lv = L * v
    
    lambdas = [10.0^k for k in 1:0.5:12]
    
    ou_vals = Float64[]
    gff_vals = Float64[]
    
    # OU Regime: A_0 = -L^{-1}
    u = L \ v
    c = dot(v, u)
    for lambda in lambdas
        Sigma_EE = Delta_E / lambda^2
        lambda_tilde = 1.0 / lambda
        L_lambda = L + lambda_tilde * I
        Lb = (c .* Lv .- v) .* Sigma_EE
        
        x1 = L_lambda \ (lambda_tilde .* Lb)
        x2 = L_lambda \ Lv
        alpha_tilde = -dot(v, x1) / dot(v, x2)
        Sigma_IE = x1 .+ alpha_tilde .* x2
        
        R_E = dot(Sigma_IE, L \ Sigma_IE)
        rho_E = R_E / Sigma_EE
        push!(ou_vals, lambda^4 * rho_E)
    end
    
    # GFF Regime: A_0 = -L
    P_I_Lv = Lv .- v .* dot(v, Lv)
    for lambda in lambdas
        Sigma_EE = Delta_E / lambda^2
        L_lambda = L + lambda * I
        
        b = P_I_Lv .* Sigma_EE
        x1 = L_lambda \ b
        x2 = L_lambda \ v
        alpha = -dot(v, x1) / dot(v, x2)
        Sigma_IE = x1 .+ alpha .* x2
        
        R_E = dot(Sigma_IE, L * Sigma_IE)
        rho_E = R_E / Sigma_EE
        push!(gff_vals, lambda^4 * rho_E)
    end
    
    # Plotting
    p = plot(lambdas, gff_vals, xscale=:log10, yscale=:log10, 
             label="GFF Regime (A = -L)", 
             xlabel="Regularization \\lambda", 
             ylabel="Trace Anomaly \\lambda^4 tr(\\rho_E)",
             title="Conditioning Spike: GFF vs OU Regime",
             linewidth=2, marker=:circle, markersize=3, color=:red,
             legend=:topleft)
             
    plot!(p, lambdas, ou_vals, 
          label="OU Regime (A = -L^{-1})", 
          linewidth=2, marker=:square, markersize=3, color=:blue)
          
    # Save figure
    mkpath("figures")
    savefig(p, "figures/gff_vs_ou_condition.png")
    println("Plot saved to figures/gff_vs_ou_condition.png")
end

run_experiment_and_plot()
