"""
    SNMRAEMSynthetics

Build a synthetic groundwater section from hydrogeological intervals and derive
everything the geophysics sees from it.

You describe the section as a list of units, each with a porosity, saturation,
clay fraction and pore-water conductivity. Two rock-physics relations then give
the two observables:

    θ_SMR   = φ · S · (1 − clay)^p_clay        free water seen by surface NMR
    σ_bulk  = σ_w · θ_total^m + σ_clay · clay  conductivity seen by airborne EM

Boundaries are gradational, porosity and clay carry correlated depth
variability, and an optional water table imposes residual saturation above it
with a capillary fringe beneath.
"""
module SNMRAEMSynthetics

using TOML, Random, Statistics, LinearAlgebra, Printf

export Unit, Column, load_config, build_column, save_column, unit_mean

"One hydrogeological interval. `bot` is filled in from the next unit's top."
struct Unit
    top::Float64
    bot::Float64
    name::String
    φ::Float64
    S::Float64
    clay::Float64
    σw::Float64
end

"A built synthetic column on a fine depth grid."
struct Column
    z::Vector{Float64}          # cell centres (m)
    dz::Float64
    φ::Vector{Float64}
    S::Vector{Float64}
    clay::Vector{Float64}
    σw::Vector{Float64}
    θ_total::Vector{Float64}
    θ_SMR::Vector{Float64}      # free water — what SMR sees
    σ_bulk::Vector{Float64}     # bulk conductivity — what AEM sees
    units::Vector{Unit}
    trans::Vector{Float64}      # transition thickness at each boundary (m)
    petro::NamedTuple           # the constants actually used, defaults resolved
    cfg::Dict{String,Any}
end

# ── config ──────────────────────────────────────────────────────────────────

"Read a TOML scenario file and return (units, cfg)."
function load_config(path::AbstractString)
    cfg = TOML.parsefile(path)
    haskey(cfg, "unit") || error("$path: no [[unit]] entries")
    raw = cfg["unit"]
    zmax = get(get(cfg, "grid", Dict()), "zmax", 200.0)
    tops = Float64[u["top"] for u in raw]
    issorted(tops) || error("[[unit]] entries must be in order of increasing top depth")
    units = Unit[]
    for (i, u) in enumerate(raw)
        bot = i < length(raw) ? raw[i+1]["top"] : zmax
        for k in ("phi", "clay", "sigma_w")
            haskey(u, k) || error("unit '$(get(u,"name","?"))' is missing `$k`")
        end
        push!(units, Unit(u["top"], bot, get(u, "name", "unit $i"),
                          u["phi"], get(u, "S", 1.0), u["clay"], u["sigma_w"]))
    end
    units, cfg
end

# ── depth-correlated variability ────────────────────────────────────────────

"Smooth zero-mean field with Gaussian correlation `corr` metres."
function gp_draw(z; corr = 6.0, seed = 1)
    n = length(z); Random.seed!(seed)
    K = [exp(-0.5 * ((z[i] - z[j]) / corr)^2) for i in 1:n, j in 1:n] + 1e-8I
    cholesky(Symmetric(K)).L * randn(n)
end

"""
Unit property at each depth, blended linearly across each gradational boundary.
`get` pulls the property off a Unit.
"""
function unit_field(z, units, trans, get)
    v = similar(z)
    @inbounds for (i, zi) in enumerate(z)
        k = findlast(u -> zi >= u.top, units); k = k === nothing ? 1 : k
        val = get(units[k])
        if k < length(units)                       # blend down into the next unit
            b = units[k+1].top; h = trans[k] / 2
            if zi > b - h
                f = clamp((zi - (b - h)) / (2h), 0, 1)
                val = (1 - f) * val + f * get(units[k+1])
            end
        end
        if k > 1                                   # blend up from the previous
            b = units[k].top; h = trans[k-1] / 2
            if zi < b + h
                f = clamp((zi - (b - h)) / (2h), 0, 1)
                val = (1 - f) * get(units[k-1]) + f * get(units[k])
            end
        end
        v[i] = val
    end
    v
end

"""
Saturation profile. With a `[watertable]` block, saturation is `s_residual`
through the vadose zone, rises across a capillary fringe of `fringe` metres and
is 1.0 at and below the water table. Without one, each unit's own `S` is used.
"""
function sat_field(z, units, trans, wt)
    isnothing(wt) && return unit_field(z, units, trans, u -> u.S)
    d  = wt["depth"]; fr = get(wt, "fringe", 0.8); sr = get(wt, "s_residual", 0.10)
    [zi >= d ? 1.0 :
     zi >= d - fr ? sr + (1 - sr) * (zi - (d - fr)) / fr : sr for zi in z]
end

# ── the build ───────────────────────────────────────────────────────────────

"""
    build_column(units, cfg) -> Column

Apply the gradational boundaries, the correlated variability and the two
rock-physics relations.
"""
function build_column(units::Vector{Unit}, cfg::Dict)
    g  = get(cfg, "grid", Dict()); dz = get(g, "dz", 0.25); zmax = get(g, "zmax", 200.0)
    z  = collect(dz/2 : dz : zmax)
    p  = get(cfg, "petrophysics", Dict())
    p_clay = get(p, "p_clay", 1.5); m = get(p, "m_archie", 1.7)
    σ_clay = get(p, "sigma_clay", 0.5)
    v  = get(cfg, "variability", Dict())
    seed = get(v, "seed", 1); corr = get(v, "corr_length", 6.0)
    sdφ = get(v, "log_sd_phi", 0.10); sdc = get(v, "log_sd_clay", 0.15)
    tmin = get(v, "transition_min", 3.0); tmax = get(v, "transition_max", 6.0)

    Random.seed!(seed)
    trans = tmin .+ (tmax - tmin) .* rand(max(length(units) - 1, 0))

    φ0   = unit_field(z, units, trans, u -> u.φ)
    clay0 = unit_field(z, units, trans, u -> u.clay)
    σw   = unit_field(z, units, trans, u -> u.σw)
    S    = sat_field(z, units, trans, get(cfg, "watertable", nothing))

    φ    = clamp.(φ0    .* exp.(gp_draw(z; corr = corr, seed = seed + 1) * sdφ), 0.005, 0.5)
    clay = clamp.(clay0 .* exp.(gp_draw(z; corr = corr, seed = seed + 2) * sdc), 0.0, 1.0)

    θ_total = φ .* S
    θ_SMR   = θ_total .* (1 .- clay) .^ p_clay
    σ_bulk  = σw .* θ_total .^ m .+ σ_clay .* clay

    Column(z, dz, φ, S, clay, σw, θ_total, θ_SMR, σ_bulk, units, trans,
           (p_clay = p_clay, m_archie = m, sigma_clay = σ_clay), cfg)
end

build_column(path::AbstractString) = build_column(load_config(path)...)

# ── reporting and output ────────────────────────────────────────────────────

"Mean of `x` over a unit, measured clear of the gradational boundaries."
function unit_mean(c::Column, x::Vector{Float64}, u::Unit; trim = 3.0)
    k = findall(zi -> u.top + trim <= zi <= u.bot - trim, c.z)
    isempty(k) && (k = findall(zi -> u.top <= zi <= u.bot, c.z))
    isempty(k) ? NaN : mean(x[k])
end

"Print the section as built, with the derived observables per unit."
function Base.show(io::IO, c::Column)
    name = get(c.cfg, "name", "unnamed")
    @printf(io, "synthetic column \"%s\" — %d units, %.2f m grid to %.0f m\n",
            name, length(c.units), c.dz, c.z[end])
    @printf(io, "petrophysics: p_clay %.2f, m %.2f, sigma_clay %.2f S/m\n",
            c.petro.p_clay, c.petro.m_archie, c.petro.sigma_clay)
    @printf(io, "%-34s %6s %6s %6s %8s %9s %8s\n",
            "unit", "phi", "S", "clay", "sigma_w", "theta_SMR", "sigma_b")
    for u in c.units
        u.bot - u.top < 2 && continue
        @printf(io, "%-34s %6.2f %6.2f %6.2f %8.3f %9.3f %8.4f\n", u.name,
                unit_mean(c, c.φ, u), unit_mean(c, c.S, u), unit_mean(c, c.clay, u),
                unit_mean(c, c.σw, u), unit_mean(c, c.θ_SMR, u), unit_mean(c, c.σ_bulk, u))
    end
end

"Write the column to JLD2 for the forward-modelling and inversion steps."
function save_column(c::Column, path::AbstractString)
    @eval Main using JLD2
    Base.invokelatest(Main.jldsave, path;
        z = c.z, dz = c.dz, φ = c.φ, S = c.S, clay = c.clay, σw = c.σw,
        θ_total = c.θ_total, θ_SMR = c.θ_SMR, σ_bulk = c.σ_bulk,
        units = [(u.top, u.bot, u.name) for u in c.units],
        trans = c.trans, config = c.cfg,
        p_clay = c.petro.p_clay, m_archie = c.petro.m_archie,
        sigma_clay = c.petro.sigma_clay)
    path
end

end # module
