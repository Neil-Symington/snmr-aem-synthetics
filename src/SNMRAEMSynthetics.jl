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

include("Lithology.jl")

export Unit, Column, load_config, build_column, save_column, unit_mean,
       check_column, lithology_names, mgL_to_sigma, sigma_to_mgL

"""
One hydrogeological interval.

Geometry can be given as `top` (bottom taken from the next unit) or as
`thickness` (tops accumulated down the section), or both if they agree. The
`sd_*` fields are per-unit log-normal spreads for the depth-correlated jitter —
a clean sand is far more uniform than a weathered profile, so one global number
for the whole section is usually wrong.
"""
struct Unit
    top::Float64
    bot::Float64
    name::String
    φ::Float64
    S::Float64
    clay::Float64
    σw::Float64
    sd_φ::Float64
    sd_clay::Float64
    sd_σw::Float64
end
thickness(u::Unit) = u.bot - u.top

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

"""
    load_config(path) -> (units, cfg)

Read a TOML scenario. Each `[[unit]]` may give its material as `lithology`
(filling porosity and clay from the preset table), its salinity as `sigma_w`
in S/m or `sigma_w_mgL` in mg/L, and its geometry as `top`, `thickness`, or
both. Anything stated explicitly wins over the preset.
"""
function load_config(path::AbstractString)
    cfg = TOML.parsefile(path)
    haskey(cfg, "unit") || error("$path: no [[unit]] entries")
    raw  = cfg["unit"]
    zmax = get(get(cfg, "grid", Dict()), "zmax", 200.0)
    v    = get(cfg, "variability", Dict())
    gsdφ = get(v, "log_sd_phi", 0.10)
    gsdc = get(v, "log_sd_clay", 0.15)
    gsdw = get(v, "log_sd_sigma_w", 0.0)   # off by default: keeps σ_w exactly as stated

    # ── geometry: resolve tops and bottoms from `top`, `thickness`, or both
    n = length(raw)
    tops = Vector{Float64}(undef, n); bots = Vector{Float64}(undef, n)
    cursor = get(raw[1], "top", 0.0)
    for (i, u) in enumerate(raw)
        nm = get(u, "name", "unit $i")
        t  = get(u, "top", cursor)
        if haskey(u, "top") && abs(t - cursor) > 1e-6 && i > 1
            error("unit \"$nm\": top $(t) m does not meet the previous unit's " *
                  "base at $(round(cursor, digits=3)) m — units must be contiguous")
        end
        th = get(u, "thickness", nothing)
        b  = th === nothing ? nothing : t + th
        tops[i] = t; bots[i] = b === nothing ? NaN : b
        cursor  = b === nothing ? NaN : b
        if b === nothing && i < n
            haskey(raw[i+1], "top") || error("unit \"$nm\": give it a `thickness`, " *
                "or give the next unit a `top` — otherwise its base is undefined")
            cursor = raw[i+1]["top"]; bots[i] = cursor
        end
    end
    isnan(bots[n]) && (bots[n] = zmax)          # last unit runs to the model base
    bots[n] > zmax && error("the section is $(bots[n]) m deep but grid.zmax is $zmax m")
    for i in 1:n
        bots[i] > tops[i] || error("unit \"$(get(raw[i],"name","unit $i"))\" has " *
                                   "non-positive thickness ($(tops[i]) to $(bots[i]) m)")
    end

    # ── properties
    units = Unit[]
    for (i, u) in enumerate(raw)
        nm   = get(u, "name", "unit $i")
        pre  = haskey(u, "lithology") ? lithology(u["lithology"]) :
               (φ = NaN, clay = NaN, sd_φ = gsdφ, sd_clay = gsdc)
        φ    = get(u, "phi",  pre.φ)
        clay = get(u, "clay", pre.clay)
        isnan(φ)    && error("unit \"$nm\": needs `phi` or a `lithology`")
        isnan(clay) && error("unit \"$nm\": needs `clay` or a `lithology`")
        σw = haskey(u, "sigma_w")     ? u["sigma_w"] :
             haskey(u, "sigma_w_mgL") ? mgL_to_sigma(u["sigma_w_mgL"]) :
             error("unit \"$nm\": needs `sigma_w` (S/m) or `sigma_w_mgL`")
        haskey(u, "sigma_w") && haskey(u, "sigma_w_mgL") &&
            error("unit \"$nm\": give `sigma_w` or `sigma_w_mgL`, not both")
        0 < φ < 1     || error("unit \"$nm\": phi $φ is not in (0,1)")
        0 <= clay <= 1 || error("unit \"$nm\": clay $clay is not in [0,1]")
        σw > 0        || error("unit \"$nm\": sigma_w must be positive")
        push!(units, Unit(tops[i], bots[i], nm, φ, get(u, "S", 1.0), clay, σw,
                          get(u, "sd_phi",     haskey(u,"lithology") ? pre.sd_φ    : gsdφ),
                          get(u, "sd_clay",    haskey(u,"lithology") ? pre.sd_clay : gsdc),
                          get(u, "sd_sigma_w", gsdw)))
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
    tmin = get(v, "transition_min", 3.0); tmax = get(v, "transition_max", 6.0)

    Random.seed!(seed)
    trans = tmin .+ (tmax - tmin) .* rand(max(length(units) - 1, 0))

    φ0   = unit_field(z, units, trans, u -> u.φ)
    clay0 = unit_field(z, units, trans, u -> u.clay)
    σw   = unit_field(z, units, trans, u -> u.σw)
    S    = sat_field(z, units, trans, get(cfg, "watertable", nothing))

    # Per-unit spreads, carried down the section as a depth profile, so a clean
    # sand can be uniform while a weathered profile is not. One correlated field
    # per property, scaled by the local spread.
    sdφ_z = unit_field(z, units, trans, u -> u.sd_φ)
    sdc_z = unit_field(z, units, trans, u -> u.sd_clay)
    sdw_z = unit_field(z, units, trans, u -> u.sd_σw)

    φ    = clamp.(φ0    .* exp.(gp_draw(z; corr = corr, seed = seed + 1) .* sdφ_z), 0.005, 0.5)
    clay = clamp.(clay0 .* exp.(gp_draw(z; corr = corr, seed = seed + 2) .* sdc_z), 0.0, 1.0)
    σw   = σw           .* exp.(gp_draw(z; corr = corr, seed = seed + 3) .* sdw_z)

    θ_total = φ .* S
    θ_SMR   = θ_total .* (1 .- clay) .^ p_clay      # free water: NMR sees this
    σ_bulk  = σw .* θ_total .^ m .+ σ_clay .* clay  # bulk conductivity: AEM sees this

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


# ── is this section a fair test? ────────────────────────────────────────────

"""
    check_column(c; theta_bar = 0.15, io = stdout) -> Vector{String}

Report whether the section is a meaningful test and flag the mistakes that look
fine on paper. Returns the list of warnings so a script can fail on them.

What it reports per unit: the split of bulk conductivity into the part carried
by pore water and the part carried by clay surfaces — which is what decides
whether salinity is recoverable at all — and whether the unit clears the
extractable bar on free water.
"""
function check_column(c::Column; theta_bar = 0.15, io = stdout)
    warn = String[]
    σc, m = c.petro.sigma_clay, c.petro.m_archie
    us = [u for u in c.units if thickness(u) >= 2.0]

    println(io, "\nconductivity budget and extractability")
    @printf(io, "%-30s %9s %9s %8s %9s %12s\n",
            "unit", "sigma_b", "clay term", "clay %", "theta_SMR", "extractable")
    info = NamedTuple[]
    for u in us
        σb = unit_mean(c, c.σ_bulk, u); θ = unit_mean(c, c.θ_SMR, u)
        cl = unit_mean(c, c.clay, u);   sw = unit_mean(c, c.σw, u)
        ct = σc * cl; pct = 100ct / σb
        ex = θ > theta_bar
        @printf(io, "%-30s %9.4f %9.4f %7.0f%% %9.3f %12s\n",
                u.name, σb, ct, pct, θ, ex ? "yes" : "no")
        push!(info, (u = u, σb = σb, θ = θ, clay = cl, σw = sw, claypct = pct, ex = ex))
    end

    # salinity is only recoverable where the water carries a real share of σ_b
    for r in info
        r.ex || continue
        # High clay share is necessary but not sufficient for the salinity to be
        # unrecoverable — it can still work if the conductivity itself is well
        # resolved. State the sensitivity rather than predict the outcome.
        amp = 1 / unit_mean(c, c.φ, r.u)^m
        r.claypct > 75 && push!(warn,
            "$(r.u.name): clay carries $(round(Int, r.claypct)) % of the conductivity, " *
            "leaving $(round(100 - r.claypct, digits = 1)) % to carry the salinity. " *
            "Recovering sigma_w divides that residual by phi^m, amplifying any error " *
            "$(round(amp, digits = 1))x, so the answer will be sensitive both to " *
            "conductivity error and to the assumed clay mineralogy (see examples/03)")
    end

    # adjacent extractable units with nothing between them are ONE aquifer
    for i in 1:length(info)-1
        info[i].ex && info[i+1].ex && push!(warn,
            "$(info[i].u.name) and $(info[i+1].u.name) both clear the " *
            "theta > $theta_bar bar and are adjacent — they will be detected as a " *
            "single unit, not two. Put a low-porosity or clay-rich interval between them.")
    end

    # a salinity step across permeable ground with no aquitard is not hydrogeology
    for i in 1:length(info)-1
        a, b = info[i], info[i+1]
        (a.ex && b.ex) || continue
        ratio = max(a.σw, b.σw) / min(a.σw, b.σw)
        ratio > 3 && push!(warn,
            "pore-water conductivity steps $(round(ratio, digits = 1))x between " *
            "$(a.u.name) and $(b.u.name), both permeable and in contact — " *
            "a step that large normally needs a confining layer to maintain")
    end

    # salinity usually increases down a flow path
    for i in 1:length(info)-1
        info[i].σw > info[i+1].σw + 1e-9 && push!(warn,
            "pore water freshens with depth from $(info[i].u.name) to " *
            "$(info[i+1].u.name) ($(round(info[i].σw, digits=3)) to " *
            "$(round(info[i+1].σw, digits=3)) S/m) — possible, but check it is intended")
    end

    # is the target enclosed by conductors? that is what decides water quality
    for i in 2:length(info)-1
        r = info[i]; r.ex || continue
        above, below = info[i-1].σb, info[i+1].σb
        (above > 2r.σb && below > 2r.σb) && push!(warn,
            "$(r.u.name) is a resistor enclosed by conductors " *
            "($(round(above, digits=3)) above, $(round(below, digits=3)) below, " *
            "unit itself $(round(r.σb, digits=3)) S/m) — the hard case: airborne EM " *
            "resolves such a unit poorly and its salinity may be set by the " *
            "regularisation rather than the data")
    end

    if isempty(warn)
        println(io, "\nno warnings.")
    else
        println(io, "\n", length(warn), " warning", length(warn) == 1 ? "" : "s", ":")
        for w in warn; println(io, "  * ", w); end
    end
    warn
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
