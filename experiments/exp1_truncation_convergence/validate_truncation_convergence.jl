"""
Experiment 1 (B.1): Truncation convergence validation.

The experiment evaluates a triangular trace-class OU witness. Truncation
diagnostics use Float64 dense linear algebra at the original N_ref = 1500
reference size.
"""

using CairoMakie
using CSV
using DataFrames
using LaTeXStrings
using LinearAlgebra
using Printf
using SHA

const OUTPUT_DIR = joinpath(@__DIR__, "results")
const NS = [16, 32, 64, 128, 256, 384, 512, 768, 1000]
const LAMBDA_GRID = Float64[2^power for power in 1:12]
const FIT_LAMBDAS = LAMBDA_GRID[end-4:end]
const SCALING_NS = [8, 32, 128, 384, 1500]

struct Model{T<:AbstractFloat}
    n_ref::Int
    gap::T
    nu::T
    coupling::T
    coupling_decay::T
    noise_decay::T
    epsilon0::T
    lambda_fixed::T
end

function model(::Type{T}; n_ref::Int = 1500) where {T<:AbstractFloat}
    return Model{T}(
        n_ref,
        T(1),
        parse(T, "0.01"),
        T(12),
        T(3),
        parse(T, "1.5"),
        T(1),
        T(4),
    )
end

function model_vectors(n::Int, m::Model{T}) where {T<:AbstractFloat}
    modes = T.(1:n)
    rates = m.gap .+ m.nu .* modes .^ 2
    outgoing = zeros(T, n)
    outgoing[2:end] .= m.coupling .* modes[2:end] .^ (-m.coupling_decay)
    noise = modes .^ (-m.noise_decay)
    return rates, outgoing, noise
end

function stationary_covariance(
    n::Int,
    m::Model{T};
    intervention_lambda::Union{Nothing,T} = nothing,
) where {T<:AbstractFloat}
    rates, outgoing, noise = model_vectors(n, m)
    if intervention_lambda !== nothing
        rates[1] = intervention_lambda
        noise[1] = m.epsilon0 / intervention_lambda
    end

    sigma = zeros(T, n, n)
    sigma[1, 1] = noise[1] / (T(2) * rates[1])
    n == 1 && return sigma

    free_rates = @view rates[2:end]
    free_outgoing = @view outgoing[2:end]
    cross = free_outgoing .* sigma[1, 1] ./ (free_rates .+ rates[1])
    sigma[2:end, 1] .= cross
    sigma[1, 2:end] .= cross

    forcing = free_outgoing * transpose(cross) + cross * transpose(free_outgoing)
    for index in eachindex(free_rates)
        forcing[index, index] += noise[index + 1]
    end
    sigma[2:end, 2:end] .= forcing ./ (free_rates .+ transpose(free_rates))
    return sigma
end

function lyapunov_max_residual(
    sigma::Matrix{T},
    m::Model{T};
    intervention_lambda::Union{Nothing,T} = nothing,
) where {T<:AbstractFloat}
    rates, outgoing, noise = model_vectors(size(sigma, 1), m)
    if intervention_lambda !== nothing
        rates[1] = intervention_lambda
        noise[1] = m.epsilon0 / intervention_lambda
    end
    residual = -(rates .+ transpose(rates)) .* sigma
    residual .+= outgoing * transpose(@view sigma[1, :])
    residual .+= (@view sigma[:, 1]) * transpose(outgoing)
    for index in eachindex(noise)
        residual[index, index] += noise[index]
    end
    return maximum(abs, residual)
end

function schur_scalar(sigma::Matrix{T}) where {T<:AbstractFloat}
    size(sigma, 1) == 1 && return sigma[1, 1]
    cross = Vector(@view sigma[2:end, 1])
    free = Symmetric(Matrix(@view sigma[2:end, 2:end]))
    return sigma[1, 1] - dot(cross, free \ cross)
end

function schur_remainder(sigma::Matrix{T}) where {T<:AbstractFloat}
    size(sigma, 1) == 1 && return zero(T)
    cross = Vector(@view sigma[2:end, 1])
    free = Symmetric(Matrix(@view sigma[2:end, 2:end]))
    return dot(cross, free \ cross)
end

function embedded_trace_distance(small::Matrix{Float64}, reference::Matrix{Float64})
    delta = -copy(reference)
    n = size(small, 1)
    delta[1:n, 1:n] .+= small
    return sum(abs, eigvals(Symmetric(delta)))
end

# Euler-Maclaurin evaluation of sum_{j=n+1}^infinity j^(-p).
function hurwitz_tail(p::Float64, n::Int)
    bernoulli = Float64[1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730]
    cutoff_terms = 32
    a = Float64(n + 1)
    x = a + cutoff_terms
    total = sum((a + index)^(-p) for index in 0:cutoff_terms-1)
    total += x^(1 - p) / (p - 1) + 0.5 * x^(-p)
    for r in eachindex(bernoulli)
        rising = prod(p + offset for offset in 0:(2 * r - 2))
        total += bernoulli[r] / factorial(2 * r) * rising * x^(-p - 2 * r + 1)
    end
    return total
end

function fitted_slope(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:AbstractFloat}
    logx = log.(x)
    logy = log.(y)
    centered_x = logx .- sum(logx) / length(logx)
    centered_y = logy .- sum(logy) / length(logy)
    return dot(centered_x, centered_y) / dot(centered_x, centered_x)
end

function intercept_in_lambda_minus_four(
    n::Int,
    m::Model{T},
    lambdas::Vector{T},
) where {T<:AbstractFloat}
    x = lambdas .^ (-4)
    y = T[
        lambda^2 * schur_scalar(
            stationary_covariance(n, m; intervention_lambda = lambda)
        ) for lambda in lambdas
    ]
    xbar = sum(x) / length(x)
    ybar = sum(y) / length(y)
    slope = dot(x .- xbar, y .- ybar) / dot(x .- xbar, x .- xbar)
    return ybar - slope * xbar
end

function truncation_diagnostics(m::Model{Float64})
    sigma_ref = stationary_covariance(m.n_ref, m)
    sigma_lambda_ref = stationary_covariance(
        m.n_ref, m; intervention_lambda = m.lambda_fixed
    )
    schur_ref = schur_scalar(sigma_lambda_ref)
    residual_ref = lyapunov_max_residual(sigma_ref, m)
    residual_lambda_ref = lyapunov_max_residual(
        sigma_lambda_ref, m; intervention_lambda = m.lambda_fixed
    )

    convergence_rows = NamedTuple[]
    for n in NS
        sigma_n = stationary_covariance(n, m)
        sigma_lambda_n = stationary_covariance(
            n, m; intervention_lambda = m.lambda_fixed
        )
        schur_n = schur_scalar(sigma_lambda_n)
        m2 = intercept_in_lambda_minus_four(n, m, FIT_LAMBDAS)
        push!(
            convergence_rows,
            (
                n = n,
                noise_trace_tail_exact = hurwitz_tail(Float64(m.noise_decay), n),
                sigma_trace_error_vs_ref = embedded_trace_distance(sigma_n, sigma_ref),
                sigma_lambda_trace_error_vs_ref = embedded_trace_distance(
                    sigma_lambda_n, sigma_lambda_ref
                ),
                schur_lambda_fixed = schur_n,
                schur_abs_error_vs_ref = abs(schur_n - schur_ref),
                m2_schur_extrapolated = m2,
                m2_schur_abs_error_theory = abs(m2 - m.epsilon0 / 2),
                zero_pole_product = 0.0,
            ),
        )
    end

    scaling_rows = NamedTuple[]
    for n in SCALING_NS
        for lambda in LAMBDA_GRID
            sigma = stationary_covariance(n, m; intervention_lambda = lambda)
            schur = schur_scalar(sigma)
            remainder = schur_remainder(sigma)
            push!(
                scaling_rows,
                (
                    n = n,
                    lambda = lambda,
                    schur = schur,
                    lambda_sq_schur = lambda^2 * schur,
                    R_lambda = remainder,
                    lambda6_R_lambda = lambda^6 * remainder,
                    m2_schur_theory = m.epsilon0 / 2,
                ),
            )
        end
    end
    return (
        DataFrame(convergence_rows),
        DataFrame(scaling_rows),
        schur_ref,
        residual_ref,
        residual_lambda_ref,
    )
end

function plot_results(convergence::DataFrame, scaling::DataFrame, m::Model{Float64}, schur_ref)
    fig = Figure(size = (1050, 740))
    ax_tail = Axis(
        fig[1, 1],
        xscale = log10,
        yscale = log10,
        title = "Trace-class truncation tails",
        xlabel = "Truncation dimension n",
        ylabel = "Error",
    )
    ax_fixed = Axis(
        fig[1, 2],
        xscale = log10,
        yscale = log10,
        title = "Soft-intervention truncation error",
        xlabel = "Truncation dimension n",
        ylabel = "Error against reference",
    )
    ax_schur = Axis(
        fig[2, 1],
        title = L"Schur residual at $\lambda=4$",
        xlabel = "Truncation dimension n",
        ylabel = L"S_{1,n,\lambda}",
    )
    ax_limit = Axis(
        fig[2, 2],
        xscale = log10,
        title = "Dimension-free scaled residual",
        xlabel = L"Intervention strength $\lambda$",
        ylabel = L"\lambda^2 S_{1,n,\lambda}",
    )

    scatterlines!(ax_tail, convergence.n, convergence.noise_trace_tail_exact; marker = :circle, label = L"\Vert D-P_nDP_n\Vert_1")
    scatterlines!(ax_tail, convergence.n, convergence.sigma_trace_error_vs_ref; marker = :rect, label = L"\Vert\Sigma_n-\Sigma_{\mathrm{ref}}\Vert_1")
    axislegend(ax_tail; position = :rt)

    scatterlines!(ax_fixed, convergence.n, convergence.sigma_lambda_trace_error_vs_ref; marker = :circle, label = L"\Vert\Sigma_{n,\lambda}-\Sigma_{\mathrm{ref},\lambda}\Vert_1")
    scatterlines!(ax_fixed, convergence.n, convergence.schur_abs_error_vs_ref; marker = :rect, label = L"|S_{n,\lambda}-S_{\mathrm{ref},\lambda}|")
    axislegend(ax_fixed; position = :rt)

    scatterlines!(ax_schur, convergence.n, convergence.schur_lambda_fixed; marker = :circle)
    hlines!(ax_schur, [schur_ref]; linestyle = :dash, color = :black)

    for n in SCALING_NS
        subset = scaling[scaling.n .== n, :]
        scatterlines!(ax_limit, subset.lambda, subset.lambda_sq_schur; marker = :circle, label = "n=$(n)")
    end
    hlines!(ax_limit, [m.epsilon0 / 2]; linestyle = :dash, color = :black)
    axislegend(ax_limit; position = :rb, nbanks = 2)

    Label(fig[0, 1:2], "Experiment C.1: Truncation Convergence and Schur Pole", fontsize = 18, font = :bold)
    save(joinpath(OUTPUT_DIR, "truncation_convergence_validation.png"), fig; px_per_unit = 2)
end

function write_summary(
    convergence::DataFrame,
    schur_ref::Float64,
    residual_ref::Float64,
    residual_lambda_ref::Float64,
)
    tail = convergence[convergence.n .>= 128, :]
    final = last(convergence)
    slopes = (
        fitted_slope(Float64.(tail.n), tail.noise_trace_tail_exact),
        fitted_slope(Float64.(tail.n), tail.sigma_trace_error_vs_ref),
        fitted_slope(Float64.(tail.n), tail.sigma_lambda_trace_error_vs_ref),
        fitted_slope(Float64.(tail.n), tail.schur_abs_error_vs_ref),
    )
    lines = [
        "# Experiment C.1 Validation and Truncation Convergence Summary",
        "",
        "## Claim Tested",
        "",
        "This Julia experiment evaluates truncation convergence for the trace-class OU witness.",
        "",
        "## Truncation Diagnostics",
        "",
        "| Quantity | Value |",
        "|---|---:|",
        @sprintf("| Reference Schur scalar at lambda=4 | %.12e |", schur_ref),
        @sprintf("| m2 extrapolation at n=1000 (Float64) | %.12e |", final.m2_schur_extrapolated),
        @sprintf("| Absolute error to epsilon0/2 at n=1000 | %.3e |", final.m2_schur_abs_error_theory),
        @sprintf("| Baseline Lyapunov max residual at N_ref | %.3e |", residual_ref),
        @sprintf("| Intervened Lyapunov max residual at N_ref | %.3e |", residual_lambda_ref),
        "",
        "| Diagnostic slope over n >= 128 | Value |",
        "|---|---:|",
        @sprintf("| Noise trace tail | %.3f |", slopes[1]),
        @sprintf("| Baseline covariance trace error | %.3f |", slopes[2]),
        @sprintf("| Intervened covariance trace error | %.3f |", slopes[3]),
        @sprintf("| Schur scalar error | %.3f |", slopes[4]),
        "",
    ]
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
    m = model(Float64)
    convergence, scaling, schur_ref, residual_ref, residual_lambda_ref =
        truncation_diagnostics(m)
    CSV.write(joinpath(OUTPUT_DIR, "truncation_convergence.csv"), convergence)
    CSV.write(joinpath(OUTPUT_DIR, "lambda_schur_scaling.csv"), scaling)
    plot_results(convergence, scaling, m, schur_ref)
    write_summary(
        convergence,
        schur_ref,
        residual_ref,
        residual_lambda_ref,
    )
    append_checksums(
        joinpath(OUTPUT_DIR, "validation_summary.md"),
        [
            joinpath(OUTPUT_DIR, "truncation_convergence.csv"),
            joinpath(OUTPUT_DIR, "lambda_schur_scaling.csv"),
            joinpath(OUTPUT_DIR, "truncation_convergence_validation.png")
        ]
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
