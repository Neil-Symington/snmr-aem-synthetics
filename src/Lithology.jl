# Lithology presets — indicative porosity and clay fraction by material, so a
# scenario can be written in the terms a hydrogeologist actually knows.
#
# These are starting points, not measurements. Override either value per unit
# wherever you know better; the preset only fills what you leave out.

const LITHOLOGY = Dict{String,NamedTuple{(:φ, :clay, :sd_φ, :sd_clay, :mineral),
                                         Tuple{Float64,Float64,Float64,Float64,String}}}(
    # name                        φ     clay   sd_φ  sd_clay  default clay mineral
    "gravel"            => (φ = 0.28, clay = 0.02, sd_φ = 0.08, sd_clay = 0.20, mineral = "kaolinite"),
    "clean sand"        => (φ = 0.32, clay = 0.04, sd_φ = 0.07, sd_clay = 0.20, mineral = "kaolinite"),
    "sand"              => (φ = 0.30, clay = 0.08, sd_φ = 0.10, sd_clay = 0.15, mineral = "kaolinite"),
    "silty sand"        => (φ = 0.28, clay = 0.18, sd_φ = 0.12, sd_clay = 0.15, mineral = "illite"),
    "sandy clay"        => (φ = 0.35, clay = 0.40, sd_φ = 0.12, sd_clay = 0.12, mineral = "illite"),
    "clay"              => (φ = 0.42, clay = 0.60, sd_φ = 0.10, sd_clay = 0.10, mineral = "illite"),
    "weathered bedrock" => (φ = 0.12, clay = 0.30, sd_φ = 0.25, sd_clay = 0.25, mineral = "kaolinite"),
    "saprock"           => (φ = 0.06, clay = 0.12, sd_φ = 0.30, sd_clay = 0.30, mineral = "kaolinite"),
    "fresh bedrock"     => (φ = 0.02, clay = 0.03, sd_φ = 0.30, sd_clay = 0.30, mineral = "kaolinite"),
)

# ── clay mineralogy ──────────────────────────────────────────────────────────
# Clay conducts at its grain surfaces, and how well depends on cation exchange
# capacity, which spans more than an order of magnitude across the common clay
# minerals. A single surface conductivity for a whole section therefore gives a
# clean kaolinitic sand the same surface conduction as a smectitic lacustrine
# clay, which is wrong by about 10x.
#
# CEC in meq/100 g, with the surface conductivity each implies:
const MINERALOGY = Dict{String,NamedTuple{(:cec_lo, :cec_hi, :σ_clay),NTuple{3,Float64}}}(
    "kaolinite" => (cec_lo =  3.0, cec_hi =  15.0, σ_clay = 0.15),
    "illite"    => (cec_lo = 10.0, cec_hi =  40.0, σ_clay = 0.50),
    "smectite"  => (cec_lo = 80.0, cec_hi = 150.0, σ_clay = 2.00),
    "mixed"     => (cec_lo = 10.0, cec_hi =  80.0, σ_clay = 0.80),
)

mineral_names() = sort(collect(keys(MINERALOGY)))

"""
    sigma_clay_for(mineral)

Surface conductivity (S/m per unit clay fraction) implied by a clay mineral.
These are indicative values tied to typical CEC ranges, NOT measurements — CEC
is rarely measured, and natural gamma (the usual clay log) responds to illite
far more than to kaolinite or smectite, so clay TYPE is normally inferred from
depositional environment rather than known. Set `sigma_clay` on a unit directly
wherever you have something better.
"""
function sigma_clay_for(mineral::AbstractString)
    key = lowercase(strip(mineral))
    haskey(MINERALOGY, key) || error("unknown clay mineral \"$mineral\". Available: " *
                                     join(mineral_names(), ", "))
    MINERALOGY[key].σ_clay
end

"Available lithology names, for error messages and `--list-lithologies`."
lithology_names() = sort(collect(keys(LITHOLOGY)))

function lithology(name::AbstractString)
    key = lowercase(strip(name))
    haskey(LITHOLOGY, key) || error("unknown lithology \"$name\". Available: " *
                                    join(lithology_names(), ", "))
    LITHOLOGY[key]
end

# Total dissolved solids to fluid conductivity. The 0.64 factor is the usual
# approximation (TDS mg/L ≈ 0.64 × EC µS/cm) and varies with ionic composition;
# good to perhaps ±15 % for ordinary groundwater.
const MGL_PER_SM = 6400.0
mgL_to_sigma(mgL) = mgL / MGL_PER_SM
sigma_to_mgL(σ)   = σ * MGL_PER_SM
