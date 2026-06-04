"""
Experiment 7 (C.3): Gaussian finite-part functional in a three-node chain.

The experiment evaluates a three-node Gaussian SDE: X_1 -> X_2 -> X_3
under soft causal intervention on X_2 with clamp strength lambda.
We compute the Schur residual and the unified renormalized clamped action J_lambda
under three cases:
  - Case A: Decoupled-centered (a_32 = 0, b_2 = 0) => J_lambda = 0 exact
  - Case B: Coupled-centered (a_32 = 1, b_2 = 0) => J_lambda = O(lambda^-4) -> 0
  - Case C: Coupled-off-centered (a_32 = 1, b_2 = 1) => J_lambda -> 1 = b_2^2/eps0
"""

using CairoMakie
using LaTeXStrings
using DataFrames
using CSV
using Printf
using SHA

const OUTPUT_DIR = joinpath(@__DIR__, "results")
const PRECISION_BITS = 256
const PLOT_POWERS = 1:12  # lambda in 2^1 to 2^12

function evaluate_chain(lam::BigFloat, a_32::BigFloat, b_2::BigFloat, eps0::BigFloat)
    # Drift matrix A_lambda:
    # A_lambda = [-1 0 0; 0 -lam 0; 0 a_32 -1]
    # Diffusion matrix D_lambda:
    # D_lambda = [1 0 0; 0 eps0/lam 0; 0 0 1]
    
    # Algebraic covariance solutions:
    sigma_11 = BigFloat(0.5)
    sigma_12 = BigFloat(0.0)
    sigma_13 = BigFloat(0.0)
    sigma_22 = eps0 / (BigFloat(2) * lam^2)
    sigma_23 = a_32 * sigma_22 / (lam + 1)
    sigma_33 = BigFloat(0.5) + a_32 * sigma_23
    
    # Schur residual S_{2, lambda}
    R_lambda = sigma_23^2 / sigma_33
    S_2_lambda = sigma_22 - R_lambda
    rho_lambda = R_lambda / sigma_22
    
    # Mean vector elements:
    m_lambda_2 = b_2 / lam
    m_lambda_3 = a_32 * m_lambda_2
    
    # Conditional bulk mean:
    m_obs_3 = a_32 * b_2 / (lam + 1)
    
    # Conditional bulk covariance:
    sigma_33_cond = sigma_33 * (1 - rho_lambda)
    
    # Bulk relative entropy:
    # Since X_1 is completely independent and identical, KL = KL_3
    # log(0) is handled
    kl_bulk = if rho_lambda == 0
        BigFloat(0.0)
    else
        BigFloat(0.5) * (sigma_33_cond / sigma_33 - 1 - log(sigma_33_cond / sigma_33) + (m_obs_3 - m_lambda_3)^2 / sigma_33)
    end
    
    # Conditional action term:
    cond_action = if S_2_lambda == 0
        BigFloat(0.0)
    else
        BigFloat(0.5) / S_2_lambda * (m_lambda_2^2 * (1 - rho_lambda)^2 + R_lambda * (1 - rho_lambda))
    end
    
    J_lambda = kl_bulk + cond_action
    
    return (
        lambda = lam,
        sigma_22 = sigma_22,
        S_2_lambda = S_2_lambda,
        lambda_sq_S = lam^2 * S_2_lambda,
        kl_bulk = kl_bulk,
        cond_action = cond_action,
        J_lambda = J_lambda,
    )
end

function plot_results(df::DataFrame)
    fig = Figure(size=(850, 360))
    
    # Left Axis (Panel A): Schur residual scaling
    ax1 = Axis(fig[1, 1],
        xscale = log10,
        title = "Schur residual scaling (Panel A)",
        xlabel = L"Clamp strength $\lambda$",
        ylabel = L"$\lambda^2 S_{2,\lambda}$",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1)
    )
    
    # Right Axis (Panel B): Renormalized action
    ax2 = Axis(fig[1, 2],
        xscale = log10,
        yscale = log10,
        title = "Renormalized clamped action (Panel B)",
        xlabel = L"Clamp strength $\lambda$",
        ylabel = L"\mathfrak{J}_\lambda",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1)
    )

    cases = [
        ("A", "Case A: Decoupled-centered", "#0B7285", :circle),
        ("B", "Case B: Coupled-centered", "#C92A2A", :rect),
        ("C", "Case C: Coupled-off-centered", "#E67700", :utriangle)
    ]

    for (c, label, color, marker) in cases
        subset = df[df.case .== c, :]
        lambda_plot = Float64.(subset.lambda)
        scaled_S = Float64.(subset.lambda_sq_S)
        # Clamp J_lambda to 1e-20 for log plot
        J_plot = [Float64(max(j, BigFloat("1e-20"))) for j in subset.J_lambda]
        
        # Plot Panel A
        lines!(ax1, lambda_plot, scaled_S, color = color, linewidth = 1.8, label = label)
        scatter!(ax1, lambda_plot, scaled_S, color = color, marker = marker, markersize = 8)
        
        # Plot Panel B
        lines!(ax2, lambda_plot, J_plot, color = color, linewidth = 1.8, label = label)
        scatter!(ax2, lambda_plot, J_plot, color = color, marker = marker, markersize = 8)
    end
    
    # Horizontal reference lines
    hlines!(ax1, [0.5], color = "#212529", linestyle = :dash, linewidth = 1.0)
    hlines!(ax2, [1.0], color = "#212529", linestyle = :dash, linewidth = 1.0)
    
    # Add legend to ax2
    axislegend(ax2, position = :rc, framevisible = false, labelsize = 10)
    
    # Supertitle
    Label(fig[0, 1:2], "Experiment C.3: three-node SDE Gaussian chain",
        fontsize = 13, font = :bold)
        
    save(joinpath(OUTPUT_DIR, "three_node_gaussian_chain.png"), fig, px_per_unit = 2.4)
end

function write_summary(df::DataFrame)
    lines = [
        "# Experiment C.3 validation summary",
        "",
        "The experiment evaluates the three-node Gaussian SDE X_1 -> X_2 -> X_3.",
        "It validates the target Schur residual scaling lambda^2 * S_{2, lambda} -> 1/2",
        "and the renormalized clamped action J_lambda under three cases.",
        "",
        "Core evaluation precision: `$(PRECISION_BITS)` bits.",
        "",
        "| Case | description | lambda^2 * S_2(lambda=4096) | J(lambda=4096) | J(lambda=2) |",
        "|---|---|---:|---:|---:|",
    ]
    
    for c in ["A", "B", "C"]
        subset = df[df.case .== c, :]
        subset = sort(subset, :lambda)
        start_row = first(subset)
        end_row = last(subset)
        
        desc = c == "A" ? "Decoupled-centered" : (c == "B" ? "Coupled-centered" : "Coupled-off-centered")
        
        end_S = @sprintf("%.8f", end_row.lambda_sq_S)
        end_J = @sprintf("%.12g", end_row.J_lambda)
        start_J = @sprintf("%.8f", start_row.J_lambda)
        
        push!(lines, "| $c | $desc | $end_S | $end_J | $start_J |")
    end
    
    push!(lines, "")
    push!(lines, "Results show:")
    push!(lines, "- Panel A: lambda^2 * S_{2, lambda} converges to 0.5 in all cases.")
    push!(lines, "- Panel B: J_lambda is exactly 0 under Case A, decays as O(lambda^-4) to 0 under Case B,")
    push!(lines, "  and converges to 1.0 under Case C.")
    
    write(joinpath(OUTPUT_DIR, "validation_summary.md"), join(lines, "\n") * "\n")
end

function append_checksums(summary_path, file_paths)
    lines = [
        "",
        "## Output File Checksums",
        "",
        "| File | MD5 | SHA-256 |",
        "|---|---|---|",
    ]
    for path in file_paths
        filename = basename(path)
        sha = bytes2hex(sha256(read(path)))
        md = "N/A"
        try
            md = readchomp(`md5 -q $path`)
        catch
        end
        push!(lines, "| $filename | $md | $sha |")
    end
    open(summary_path, "a") do io
        write(io, join(lines, "\n") * "\n")
    end
end

function main()
    mkpath(OUTPUT_DIR)
    setprecision(BigFloat, PRECISION_BITS) do
        eps0 = BigFloat(1)
        lambdas = [BigFloat(2)^power for power in PLOT_POWERS]
        
        records = []
        # Case A: a_32 = 0, b_2 = 0
        for lam in lambdas
            rec = evaluate_chain(lam, BigFloat(0), BigFloat(0), eps0)
            push!(records, merge((case="A",), rec))
        end
        # Case B: a_32 = 1, b_2 = 0
        for lam in lambdas
            rec = evaluate_chain(lam, BigFloat(1), BigFloat(0), eps0)
            push!(records, merge((case="B",), rec))
        end
        # Case C: a_32 = 1, b_2 = 1
        for lam in lambdas
            rec = evaluate_chain(lam, BigFloat(1), BigFloat(1), eps0)
            push!(records, merge((case="C",), rec))
        end
        
        df = DataFrame(records)
        CSV.write(joinpath(OUTPUT_DIR, "three_node_data.csv"), df)
        plot_results(df)
        write_summary(df)
        append_checksums(
            joinpath(OUTPUT_DIR, "validation_summary.md"),
            [
                joinpath(OUTPUT_DIR, "three_node_data.csv"),
                joinpath(OUTPUT_DIR, "three_node_gaussian_chain.png")
            ]
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
