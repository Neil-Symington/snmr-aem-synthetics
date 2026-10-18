#!/usr/bin/env julia
# Example 3 — why conductivity alone cannot give you salinity.
#
#   julia --project=. examples/03_why_clay_matters.jl
#
# Takes the aquifer from the minimal section and splits its conductivity into
# the part carried by the pore water and the part carried by clay surfaces.
# Then asks: if you measured that conductivity perfectly but assumed the wrong
# clay mineralogy, what salinity would you infer?

include(joinpath(@__DIR__, "..", "src", "SNMRAEMSynthetics.jl"))
using .SNMRAEMSynthetics, Statistics, Printf

units, cfg = load_config(joinpath(@__DIR__, "..", "config", "minimal.toml"))
col = build_column(units, cfg)

aq  = units[1]
k   = findall(z -> aq.top + 3 <= z <= aq.bot - 3, col.z)
φ   = mean(col.φ[k]); clay = mean(col.clay[k])
σw  = mean(col.σw[k]); σb = mean(col.σ_bulk[k]); θ = mean(col.θ_total[k])
m   = col.petro.m_archie      # resolved by the builder, defaults and all
σc  = col.petro.sigma_clay

# the builder keeps the two contributions separately (works for both the
# Simandoux and the Waxman-Smits form)
water = mean(col.σ_water[k])
clayt = mean(col.σ_clay[k])

@printf("""
Aquifer unit: phi %.2f, clay %.3f, true pore water %.3f S/m

  bulk conductivity          %.4f S/m
    carried by the water     %.4f S/m   (%.0f%%)
    carried by clay surfaces %.4f S/m   (%.0f%%)

Everything you know about salinity lives in that water term. Recovering it
means subtracting the clay term and dividing by phi^m = %.4f, which multiplies
any error by %.1f.

""", φ, clay, σw, σb, water, 100water/σb, clayt, 100clayt/σb, φ^m, 1/φ^m)

println("Now assume a clay type and infer the salinity from a PERFECT conductivity:\n")
@printf("%-26s %12s %14s %12s\n", "assumed clay", "sigma_clay", "inferred sw", "x truth")
for (name, sc) in (("kaolinite (low CEC)", 0.15), ("illite", 0.30),
                   ("illite / mixed  <-- truth", 0.50), ("smectite-rich", 1.00),
                   ("smectite (high CEC)", 2.00))
    inferred = (σb - sc*clay) / φ^m
    @printf("%-26s %12.2f %14s %12s\n", name, sc,
            inferred <= 0 ? "inadmissible" : @sprintf("%.3f", inferred),
            inferred <= 0 ? "—" : @sprintf("%.1f", inferred/σw))
end

# The same question under Waxman-Smits, where the clay term is carried through
# the pore network and parameterised by CEC rather than a bare conductivity
cfg_ws = deepcopy(cfg); cfg_ws["petrophysics"] = Dict("model" => "waxman_smits", "cec_clay" => 25.0)
col_ws = build_column(units, cfg_ws)
@printf("\nWaxman-Smits (CEC of the clay fraction 25 meq/100 g, illite-like): bulk %.4f S/m, water %.0f%%, clay %.0f%%\n",
        mean(col_ws.σ_bulk[k]), 100mean(col_ws.σ_water[k])/mean(col_ws.σ_bulk[k]),
        100mean(col_ws.σ_clay[k])/mean(col_ws.σ_bulk[k]))
@printf("%-26s %12s %14s %12s\n", "assumed clay", "CEC meq/100g", "inferred sw", "x truth")
for (name, cec) in (("kaolinite", 8.0), ("illite  <-- truth", 25.0), ("mixed-layer", 50.0), ("smectite", 100.0))
    φw = mean(col_ws.φ[k]); Q = SNMRAEMSynthetics.qv(cec * clay, φw, col_ws.petro.rho_grain)
    # solve σ_b = φ^m (σ_w + B(σ_w) Q) for σ_w by fixed-point iteration
    sw = σw
    for _ in 1:50
        sw = mean(col_ws.σ_bulk[k]) / φw^m - SNMRAEMSynthetics.waxman_B(max(sw, 1e-4)) * Q
    end
    @printf("%-26s %12.0f %14s %12s\n", name, cec,
            sw <= 0 ? "inadmissible" : @sprintf("%.3f", sw), sw <= 0 ? "—" : @sprintf("%.1f", sw/σw))
end

println("""

The conductivity was exact in every row. Only the assumed mineralogy changed,
and the inferred salinity moves by nearly an order of magnitude — far more than
the fresh/brackish boundary at about 0.15 S/m (~1,000 mg/L).

This is why the synthetics exist: to find out when the geophysics can answer a
water-quality question and when it is really answering a question about your
assumptions. A gamma log gives clay CONTENT; spectral gamma or XRD is what
constrains clay TYPE.
""")
