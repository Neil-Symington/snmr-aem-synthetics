#!/usr/bin/env julia
# Forward-model the airborne EM and surface NMR responses of a built column.
#
#   julia --project=. scripts/forward.jl output/unconfined_aquifer_truth.jld2
#
# Writes <name>_forward.jld2 with the noisy AEM operator state and the SMR
# sounding curve, ready for inversion.
using JLD2, Statistics, Printf, Random, LinearAlgebra
using HiQGA.transD_GP, SNMRForward
include(joinpath(@__DIR__, "..", "src", "aem_setup.jl"))

length(ARGS) >= 1 || error("usage: forward.jl <truth.jld2>")
T = load(ARGS[1]); cfg = T["config"]
out = replace(ARGS[1], "_truth.jld2" => "_forward.jld2")

acfg = get(cfg, "aem", Dict()); scfg = get(cfg, "smr", Dict())

# ── AEM: SkyTEM-class system, AusAEM-standard noise
grid  = build_aem_grid()
σ_aem = piecewise_interp(grid.z, T["z"], T["σ_bulk"]); σ_aem[grid.z .< 0] .= 1e-12
built = build_skytem_operator(T["σ_bulk"], T["z"];
            noise_scale = get(acfg, "noise_scale", 1.0),
            rseed_noise = get(acfg, "rseed", 123))
aem = built.aem
@printf("AEM: %d LM + %d HM gates over %d cells\n",
        length(aem.dlow), length(aem.dhigh), length(grid.zall))

# ── SMR: square loop, kernel built from the TRUE conductivity
#    (inverting with a kernel from an estimated conductivity is the realistic
#     case and is what the inversion scripts do; this is the reference forward)
loop = get(scfg, "loop_side", 88.6)
Be   = get(scfg, "Be_nT", 58000.0) * 1e-9
inc  = deg2rad(get(scfg, "inclination_deg", -66.5))
zsmr = collect(get(scfg, "z_first", 1.0) : get(scfg, "dz", 2.0) : get(scfg, "z_last", 120.0))
qs   = get(scfg, "q", collect(range(0.5, 12.0, length = 24)))
# The kernel's conductivity model is carried on the AEM grid's own layering —
# interpolating onto the coarser SMR grid first throws away resolution the
# kernel can use. MRSForward_square takes the loop HALF-side.
zA   = grid.zall[grid.zall .> 0]
σA   = piecewise_interp(zA, T["z"], T["σ_bulk"])
ct   = SNMRForward.ConductivityModel(σA, diff(zA))
F    = SNMRForward.MRSForward_square(loop/2, zsmr, qs, inc, 0.0, Be, ct)
θ_on = [mean(T["θ_SMR"][findall(zi -> abs(zi - zs) <= 1.5, T["z"])]) for zs in zsmr]
# SNMRForward.forward does trapezoidal integration over depth with the volume
# element already folded into the stored kernel — it is NOT kernel * w.
V0   = abs.(SNMRForward.forward(F, θ_on))
@printf("SMR: %d pulse moments, V0 %.1f–%.1f nV\n", length(qs), 1e9*minimum(V0), 1e9*maximum(V0))

Random.seed!(get(scfg, "rseed", 47))
σ̂ = get(scfg, "noise_nV", 4.63) * 1e-9 .* (0.6 .+ 0.8 .* rand(length(qs)))
Vn = max.(V0 .+ σ̂ .* randn(length(qs)), 0.0)
@printf("SMR noise: median %.2f nV, SNR %.0f–%.0f\n",
        1e9*median(σ̂), minimum(V0 ./ σ̂), maximum(V0 ./ σ̂))

jldsave(out; zall_aem = grid.zall, σ_aem_truth = built.σ_aem_truth,
        dlow = aem.dlow, dhigh = aem.dhigh, σlow = aem.σlow, σhigh = aem.σhigh,
        zsmr = zsmr, q = qs, V0_clean = V0, V0_noisy = Vn, σ̂ = σ̂,
        loop_side = loop, Be = Be, inclination = inc, config = cfg)
println("wrote ", out)
