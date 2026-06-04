"""
Experiment 6 (C.2): exact-calibration ablation for target covariance.

The experiment evaluates the target-coordinate covariance calibration ratio and defect
under four candidate target-diffusion families.
"""

using CairoMakie
using LaTeXStrings
using DataFrames
using CSV
using Printf
using SHA

const OUTPUT_DIR = joinpath(@__DIR__, "results")
const PRECISION_BITS = 256
const PLOT_POWERS = 0:12
const PRECISION_PROBE_POWER = 100

struct ScalingFamily
    key::String
    label::String
    formula::LaTeXString
    color::String
    marker::Symbol
end

const FAMILIES = (
    ScalingFamily(
        "canonical",
        "Canonical",
        L"d_\lambda=\varepsilon_0/\lambda",
        "#0B7285",
        :circle,
    ),
    ScalingFamily(
        "slow",
        "Slow decay",
        L"d_\lambda=\varepsilon_0",
        "#C92A2A",
        :rect,
    ),
    ScalingFamily(
        "fast",
        "Fast decay",
        L"d_\lambda=\varepsilon_0/\lambda^2",
        "#E67700",
        :utriangle,
    ),
    ScalingFamily(
        "perturbed",
        "Asymptotic perturbation",
        L"d_\lambda=(\varepsilon_0/\lambda)(1+\lambda^{-2})",
        "#5F3DC4",
        :diamond,
    ),
)

function diffusion_variance(key::String, lam::BigFloat, eps0::BigFloat)
    if key == "canonical"
        return eps0 / lam
    elseif key == "slow"
        return eps0
    elseif key == "fast"
        return eps0 / lam^2
    elseif key == "perturbed"
        return (eps0 / lam) * (one(BigFloat) + lam^-2)
    else
        error("Unknown scaling family: $key")
    end
end

function evaluate_family(
    family::ScalingFamily,
    lam::BigFloat,
    eps0::BigFloat,
)
    d_lambda = diffusion_variance(family.key, lam, eps0)
    # The target Lyapunov equation gives Sigma_kk = d_lambda / (2 * lam)
    # Thus the calibration ratio 2 * lam^2 * Sigma_kk / eps0 is exactly lam * d_lambda / eps0
    ratio = lam * d_lambda / eps0
    defect = abs(ratio - 1)
    return (
        family = family.key,
        label = family.label,
        formula = string(family.formula),
        lambda = lam,
        d_lambda = d_lambda,
        ratio = ratio,
        defect = defect,
        precision_bits = PRECISION_BITS,
    )
end

function plot_results(df::DataFrame)
    # 11.2 x 4.5 inches equivalent resolution (e.g. 806 x 324 points)
    fig = Figure(size=(850, 360))
    
    # Left Axis (Panel A): Calibration ratio
    ax1 = Axis(fig[1, 1],
        xscale = log10,
        yscale = log10,
        title = "Calibration ratio (Panel A)",
        xlabel = L"Clamp strength $\lambda$",
        ylabel = L"$\lambda d_\lambda/\varepsilon_0$",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1),
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminorticks = IntervalsBetween(9),
        yminorticks = IntervalsBetween(9)
    )
    
    # Right Axis (Panel B): Exact-calibration defect
    ax2 = Axis(fig[1, 2],
        xscale = log10,
        yscale = log10,
        title = "Calibration defect (Panel B)",
        xlabel = L"Clamp strength $\lambda$",
        ylabel = L"$|\lambda d_\lambda/\varepsilon_0 - 1|$",
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:black, 0.1),
        ygridcolor = (:black, 0.1),
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminorticks = IntervalsBetween(9),
        yminorticks = IntervalsBetween(9)
    )

    for family in FAMILIES
        subset = df[df.family .== family.key, :]
        lambda_plot = Float64.(subset.lambda)
        ratio_plot = Float64.(subset.ratio)
        # Handle defect = 0 by clamping it for log plot
        defect_plot = [Float64(max(d, BigFloat("1e-20"))) for d in subset.defect]
        
        # Line + Scatter for ax1
        lines!(ax1, lambda_plot, ratio_plot,
            color = family.color, linewidth = 1.8, label = family.label)
        scatter!(ax1, lambda_plot, ratio_plot,
            color = family.color, marker = family.marker, markersize = 8)
            
        # Line + Scatter for ax2
        lines!(ax2, lambda_plot, defect_plot,
            color = family.color, linewidth = 1.8, label = family.label)
        scatter!(ax2, lambda_plot, defect_plot,
            color = family.color, marker = family.marker, markersize = 8)
    end
    
    # Horizontal reference lines
    hlines!(ax1, [1.0], color = "#212529", linestyle = :dash, linewidth = 1.0)
    
    # Add legend to ax2
    axislegend(ax2, position = :rt, framevisible = false, labelsize = 10)
    
    # Supertitle
    Label(fig[0, 1:2], "Experiment C.2: covariance calibration diagnostic",
        fontsize = 13, font = :bold)
        
    save(joinpath(OUTPUT_DIR, "ablation_combined.png"), fig, px_per_unit = 2.4)
end

function write_summary(df::DataFrame, probe_df::DataFrame)
    decimal_digits = floor(Int, PRECISION_BITS * log10(2))
    perturbed_probe = only(eachrow(probe_df[probe_df.family .== "perturbed", :]))
    probe_defect = perturbed_probe.defect
    lines = [
        "# Experiment C.2 validation summary",
        "",
        "The experiment evaluates the target covariance calibration ratio `lambda * d_lambda / eps0`",
        "and the defect `|lambda * d_lambda / eps0 - 1|`.",
        "",
        "Core evaluation precision: `$(PRECISION_BITS)` bits (approximately `$(decimal_digits)` decimal digits).",
        "",
        "| family | ratio(lambda=1) | ratio(lambda=4096) | defect(lambda=4096) | result |",
        "|---|---:|---:|---:|---|",
    ]
    results = Dict(
        "canonical" => "exactly scale invariant",
        "slow" => "ratio diverges",
        "fast" => "ratio collapses to zero",
        "perturbed" => "asymptotically correct, not scale invariant",
    )
    for family in FAMILIES
        subset = df[df.family .== family.key, :]
        subset = sort(subset, :lambda)
        start_row = first(subset)
        end_row = last(subset)
        
        start_ratio = @sprintf("%.8g", start_row.ratio)
        end_ratio = @sprintf("%.8g", end_row.ratio)
        end_defect = @sprintf("%.8g", end_row.defect)
        
        push!(lines, "| $(family.label) | $(start_ratio) | $(end_ratio) | $(end_defect) | $(results[family.key]) |")
    end
    push!(lines, "")
    push!(lines, "The asymptotic perturbation family reaches the calibrated limit but fails the")
    push!(lines, "exact finite-lambda invariance condition. The experiment checks the consequence")
    push!(lines, "of the covariance calibration condition.")
    push!(lines, "")
    push!(lines, "High-precision probe: at `lambda = 2^$(PRECISION_PROBE_POWER)`, the perturbed")
    push!(lines, "family remains non-constant with exact-calibration defect = $(probe_defect).")
    
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
        records = [
            evaluate_family(family, lam, eps0)
            for family in FAMILIES for lam in lambdas
        ]
        df = DataFrame(records)
        probe_lambda = BigFloat(2)^PRECISION_PROBE_POWER
        probe_df = DataFrame(
            evaluate_family(family, probe_lambda, eps0)
            for family in (FAMILIES[1], FAMILIES[4])
        )
        CSV.write(joinpath(OUTPUT_DIR, "ablation_data.csv"), df)
        CSV.write(joinpath(OUTPUT_DIR, "precision_probe.csv"), probe_df)
        plot_results(df)
        write_summary(df, probe_df)
        append_checksums(
            joinpath(OUTPUT_DIR, "validation_summary.md"),
            [
                joinpath(OUTPUT_DIR, "ablation_data.csv"),
                joinpath(OUTPUT_DIR, "precision_probe.csv"),
                joinpath(OUTPUT_DIR, "ablation_combined.png")
            ]
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
