#!/usr/bin/env julia
# Example 2 — the same thing without a TOML file.
#
#   julia --project=. examples/02_build_in_code.jl
#
# Sometimes you want a section generated in a loop rather than written out by
# hand — a suite of scenarios, a sensitivity sweep, a Monte Carlo over
# geometries. The config is just a Dict, so you can build it in code.

include(joinpath(@__DIR__, "..", "src", "SNMRAEMSynthetics.jl"))
using .SNMRAEMSynthetics, Printf

# Unit(top, bot, name, phi, S, clay, sigma_w; sd_phi, sd_clay, sd_sigma_w, cec)
# — the keyword arguments are optional and fall back to the section defaults
units = [
    Unit( 0.0, 25.0, "sand aquifer",  0.30, 1.0, 0.05, 0.05),
    Unit(25.0, 50.0, "clay aquitard", 0.40, 1.0, 0.60, 0.40),
    Unit(50.0, 120.0, "basement",     0.03, 1.0, 0.02, 0.80),
]

cfg = Dict{String,Any}(
    "name" => "in_code",
    "grid" => Dict("dz" => 0.25, "zmax" => 120.0),
    "petrophysics" => Dict("p_clay" => 1.5, "m_archie" => 1.7, "sigma_clay" => 0.5),
    "variability" => Dict("seed" => 42, "corr_length" => 6.0),
)

show(stdout, build_column(units, cfg))

# Vary one thing across a suite — here, how thick the aquifer is.
println("\n\nSweeping aquifer thickness:\n")
@printf("%-10s %12s %14s\n", "thickness", "theta_SMR", "sigma_bulk")
for base in (15.0, 25.0, 40.0)
    u = [Unit(0.0, base, "sand aquifer", 0.30, 1.0, 0.05, 0.05),
         Unit(base, base+25, "clay aquitard", 0.40, 1.0, 0.60, 0.40),
         Unit(base+25, 120.0, "basement", 0.03, 1.0, 0.02, 0.80)]
    c = build_column(u, cfg)
    k = findall(z -> 3 <= z <= base - 3, c.z)
    @printf("%-10s %12.3f %14.4f\n", "$(base) m",
            sum(c.θ_SMR[k])/length(k), sum(c.σ_bulk[k])/length(k))
end
