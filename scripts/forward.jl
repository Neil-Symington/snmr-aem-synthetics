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
built = build_skytem_operator(T["σ_bulk"], T["z"];
            noise_scale = get(acfg, "noise_scale", 1.0),
            noisefrac   = get(acfg, "noise_frac", 0.03),
            rseed_noise = get(acfg, "rseed", 123))
aem  = built.aem
grid = built.grid
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
# kernel can use. Layer interfaces are the AEM cell BOUNDARIES (grid.zbnd,
# which starts at 0) and each layer takes the conductivity at its cell centre.
# MRSForward_square takes the loop HALF-side (verified in SMRPInversion tests).
σA   = piecewise_interp(grid.zall, T["z"], T["σ_bulk"])
ct   = SNMRForward.ConductivityModel(σA, diff(grid.zbnd))
F    = SNMRForward.MRSForward_square(loop/2, zsmr, qs, inc, 0.0, Be, ct)
# water content is a NODAL field on zsmr (trapezoidal rule in SNMRForward.forward);
# average the fine truth over ±1.5 m around each node
θ_on = [mean(T["θ_SMR"][findall(zi -> abs(zi - zs) <= 1.5, T["z"])]) for zs in zsmr]
# SNMRForward.forward does trapezoidal integration over depth with the volume
# element already folded into the stored kernel — it is NOT kernel * w.
Vc   = SNMRForward.forward(F, θ_on)          # complex response
V0   = abs.(Vc)
@printf("SMR: %d pulse moments, V0 %.1f–%.1f nV\n", length(qs), 1e9*minimum(V0), 1e9*maximum(V0))

# Noise. The per-gate shape is random within ±40 % of noise_nV so the
# hierarchical inversion has a heteroscedastic shape to learn. Default 4.63 nV
# is a quiet site after heavy stacking; set smr.noise_nV from a measured record
# (e.g. derive it from the stack scatter of a field sounding) for realism.
#   noise_domain = "amplitude": Gaussian noise added to |V| (high-SNR limit)
#   noise_domain = "complex":   Gaussian noise on each quadrature of V, then
#                               modulus — Rician amplitudes, correct at low SNR
#                               (Rice 1944; bias +14 % at SNR 2, +55 % at SNR 1)
rng = Xoshiro(get(scfg, "rseed", 47))
σ̂ = get(scfg, "noise_nV", 4.63) * 1e-9 .* (0.6 .+ 0.8 .* rand(rng, length(qs)))
noise_domain = get(scfg, "noise_domain", "amplitude")
if noise_domain == "complex"
    Vn = abs.(Vc .+ σ̂ .* (randn(rng, length(qs)) .+ im .* randn(rng, length(qs))))
elseif noise_domain == "amplitude"
    Vn = max.(V0 .+ σ̂ .* randn(rng, length(qs)), 0.0)
else
    error("smr.noise_domain must be \"amplitude\" or \"complex\"")
end
@printf("SMR noise (%s): median %.2f nV, SNR %.0f–%.0f\n",
        noise_domain, 1e9*median(σ̂), minimum(V0 ./ σ̂), maximum(V0 ./ σ̂))

jldsave(out; zall_aem = grid.zall, zbnd_aem = grid.zbnd, σ_aem_truth = built.σ_aem_truth,
        dlow = aem.dlow, dhigh = aem.dhigh, σlow = aem.σlow, σhigh = aem.σhigh,
        zsmr = zsmr, q = qs, V0_clean = V0, V_clean_complex = Vc, V0_noisy = Vn, σ̂ = σ̂,
        noise_domain = noise_domain,
        loop_side = loop, Be = Be, inclination = inc, config = cfg)
println("wrote ", out)
