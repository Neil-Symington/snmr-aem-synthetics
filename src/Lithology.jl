# Lithology presets — indicative porosity and clay fraction by material, so a
# scenario can be written in the terms a hydrogeologist actually knows.
#
# These are starting points, not measurements. Override either value per unit
# wherever you know better; the preset only fills what you leave out.

const LITHOLOGY = Dict{String,NamedTuple{(:φ, :clay, :sd_φ, :sd_clay),NTuple{4,Float64}}}(
    # name                        φ     clay   sd_φ  sd_clay
    "gravel"            => (φ = 0.28, clay = 0.02, sd_φ = 0.08, sd_clay = 0.20),
    "clean sand"        => (φ = 0.32, clay = 0.04, sd_φ = 0.07, sd_clay = 0.20),
    "sand"              => (φ = 0.30, clay = 0.08, sd_φ = 0.10, sd_clay = 0.15),
    "silty sand"        => (φ = 0.28, clay = 0.18, sd_φ = 0.12, sd_clay = 0.15),
    "sandy clay"        => (φ = 0.35, clay = 0.40, sd_φ = 0.12, sd_clay = 0.12),
    "clay"              => (φ = 0.42, clay = 0.60, sd_φ = 0.10, sd_clay = 0.10),
    "weathered bedrock" => (φ = 0.12, clay = 0.30, sd_φ = 0.25, sd_clay = 0.25),
    "saprock"           => (φ = 0.06, clay = 0.12, sd_φ = 0.30, sd_clay = 0.30),
    "fresh bedrock"     => (φ = 0.02, clay = 0.03, sd_φ = 0.30, sd_clay = 0.30),
)

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
