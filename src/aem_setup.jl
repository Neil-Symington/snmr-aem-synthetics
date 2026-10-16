# aem_setup.jl
#
# Build a SkyTEM AEM operator for a given σ(z) truth model from a PETRO
# scenario JLD2. Uses the AusAEM-realistic electronics_halt noise pattern
# adapted from HiQGA/examples/SkyTEM1D/synth/McMC/electronics_halt.jl.
#
# Noise model (HiQGA makenoisydata!): σ_gate² = (noisefrac·|dB/dt|)² + σ_halt²,
# i.e. a multiplicative floor plus the additive late-time (electronics) noise.
# Trans-D GP prior and parallel tempering: Ray & Myer (2019), GJI 217,
# 1706–1726; Ray et al. (2013), Geophysics 78, E271–E280.

using HiQGA
using HiQGA.transD_GP, JLD2, Random, LinearAlgebra, Statistics

# Path to the electronics_halt file bundled with HiQGA — resolved via pathof(HiQGA).
# Included at top level so the constants (LM_times, HM_times, LM_ramp, HM_ramp,
# LM_noise, HM_noise, rTx, lowpassfcs) leak into scope where build_skytem_operator
# can see them.
const HIQGA_ROOT = joinpath(dirname(pathof(HiQGA)), "..")
const HIQGA_ELECTRONICS_HALT = joinpath(HIQGA_ROOT, "examples", "SkyTEM1D", "synth",
                                        "McMC", "electronics_halt.jl")
include(HIQGA_ELECTRONICS_HALT)

"""
    load_scenario_conductivity(petro_id, dir)

Load σ_bulk (S/m) and zboundaries (m) from a PETRO JLD2 file. dir is the
directory containing petro_<id>.jld2 files.
"""

"""
    piecewise_interp(z_target, z_source, v_source)

Piecewise-constant interpolation of `v_source(z_source)` onto `z_target`.
Preserves layer structure — for each target depth, uses the value at the
nearest deeper source boundary <=.
"""
function piecewise_interp(z_target::AbstractVector, z_source::AbstractVector,
                          v_source::AbstractVector)
    v = similar(z_target, Float64)
    for (i, zt) in enumerate(z_target)
        idx = searchsortedlast(z_source, zt)
        idx = clamp(idx, 1, length(v_source))
        v[i] = v_source[idx]
    end
    v
end

"""
    build_aem_grid()

AEM depth grid — log-growing from surface to ~280 m, adapted from HiQGA
SkyTEM synth example. Returns (zall, znall, zboundaries, z, ρ_air).
"""
function build_aem_grid()
    zfixed = [-1e5]
    ρfixed = [1e12]
    zstart = 0.0
    extendfrac, dz = 1.03, 1.5
    zall, znall, zbnd =
        transD_GP.setupz(zstart, extendfrac, dz=dz, n=65, showplot=false)
    z, ρ, nfixed = transD_GP.makezρ(zbnd; zfixed=zfixed, ρfixed=ρfixed)
    (zall=zall, znall=znall, zbnd=zbnd, z=z, ρ_initial=ρ, nfixed=nfixed)
end

"""
    build_skytem_operator(σ_smr, zbnd_smr; noise_scale=1.0, noisefrac=0.03, rseed_noise=47)

Build a SkyTEM `dBzdt` operator with the σ(z) truth interpolated onto the
AEM depth grid, forward-modelled, and noise added. `noise_scale` multiplies
the additive halt noise, `noisefrac` is the multiplicative floor.

Returns `(aem, aem_grid, ρ_aem_truth, σ_aem_truth)`.
"""
function build_skytem_operator(σ_smr::AbstractVector, zbnd_smr::AbstractVector;
                               noise_scale::Real = 1.0, noisefrac::Real = 0.03,
                               rseed_noise::Int = 47)
    aem_grid = build_aem_grid()

    # Interpolate truth σ onto AEM grid, invert to ρ
    σ_aem_truth = piecewise_interp(aem_grid.z, zbnd_smr, σ_smr)
    σ_aem_truth[aem_grid.z .< 0] .= 1e-12
    ρ_aem_truth = 1.0 ./ σ_aem_truth

    # SkyTEM geometry (LM_*, HM_*, rTx, lowpassfcs from top-level electronics_halt include)
    rRx = 13.0; zRx = -42.0; zTx = -40.0
    LM_noise_scaled = LM_noise .* noise_scale
    HM_noise_scaled = HM_noise .* noise_scale

    aem = transD_GP.SkyTEM1DInversion.dBzdt(;
        timeslow = LM_times, ramplow = LM_ramp, zRxlow=zRx, zTxlow = zTx,
        timeshigh = HM_times, ramphigh = HM_ramp, zRxhigh=zRx, zTxhigh = zTx,
        rRx = rRx, rTx = rTx, z = aem_grid.z, ρ = copy(ρ_aem_truth),
        lowpassfcs = lowpassfcs)

    # Forward + noise
    # makenoisydata! reseeds internally from its OWN rseed keyword (default 123),
    # so a preceding Random.seed! has no effect. It must be passed through, or
    # every "different" seed gives byte-identical data.
    transD_GP.SkyTEM1DInversion.makenoisydata!(aem, log10.(ρ_aem_truth[2:end]);
        rseed = rseed_noise, showplot = false, noisefrac = noisefrac,
        σ_halt_low = LM_noise_scaled, σ_halt_high = HM_noise_scaled, units = 1e-12)

    (aem = aem, grid = aem_grid, ρ_aem_truth = ρ_aem_truth, σ_aem_truth = σ_aem_truth)
end

"""
    make_aem_opt(znall_aem, fileprefix; nmin=2, nmax=40, fbounds=[-0.5 3.25])

Set up the transD_GP OptionsStat for AEM (log10-ρ) MCMC. Adapted from
HiQGA/examples/SkyTEM1D/synth/McMC/02_make_aem_inversion_opts_s.jl.
"""
function make_aem_opt(znall_aem::AbstractVector, fileprefix::AbstractString;
                      nmin::Int = 2, nmax::Int = 40,
                      fbounds::AbstractMatrix = [-0.5 3.25])
    K = transD_GP.GP.OrstUhn()
    sdev_pos  = [0.05*abs(diff([extrema(znall_aem)...])[1])]
    sdev_prop = 0.07*diff(fbounds, dims=2)[:]
    sdev_dc   = 0.008*diff(fbounds, dims=2)[:]
    xall      = permutedims(collect(znall_aem))
    xbounds   = permutedims([extrema(znall_aem)...])
    transD_GP.OptionsStat(
        nmin = nmin, nmax = nmax,
        xbounds = xbounds, fbounds = fbounds,
        fdataname = fileprefix, xall = xall,
        λ = [2.0], δ = 0.1,
        sdev_prop = sdev_prop, sdev_pos = sdev_pos, sdev_dc = sdev_dc,
        save_freq = 50, demean = false, sampledc = true,
        quasimultid = false, K = K,
    )
end
