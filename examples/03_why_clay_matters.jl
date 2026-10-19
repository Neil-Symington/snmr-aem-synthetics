#!/usr/bin/env julia
# Example 3 — why conductivity alone cannot give you salinity.
#
#   julia --project=. examples/03_why_clay_matters.jl
#
# Takes the aquifer from the minimal section and splits its conductivity into
# the part carried by the pore water and the part carried by clay surfaces.
# Then asks: if you measured that conductivity perfectly but assumed the wrong
# clay mineralogy, what salinity would you infer?
#
# Clay mineral is an input per unit — `clay_mineral = "kaolinite"` and so on, or
# `sigma_clay` directly. Try changing it in config/minimal.toml and re-running:
# the inferred salinity moves far more than most people expect, which is the
# whole point of the example.

include(joinpath(@__DIR__, "..", "src", "SNMRAEMSynthetics.jl"))
using .SNMRAEMSynthetics, Statistics, Printf

units, cfg = load_config(joinpath(@__DIR__, "..", "config", "minimal.toml"))
col = build_column(units, cfg)

aq  = units[1]
k   = findall(z -> aq.top + 3 <= z <= aq.bot - 3, col.z)
φ   = mean(col.φ[k]); clay = mean(col.clay[k])
σw  = mean(col.σw[k]); σb = mean(col.σ_bulk[k]); θ = mean(col.θ_total[k])
m   = col.petro.m_archie      # resolved by the builder, defaults and all
# Surface conduction is a property of THIS unit now, set by its clay mineral
# (or stated directly). It is the single input that most changes the answer.
σc  = aq.σ_clay

water = σw * θ^m
clayt = σc * clay

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

println("""

The conductivity was exact in every row. Only the assumed mineralogy changed,
and the inferred salinity moves by nearly an order of magnitude — far more than
the fresh/brackish boundary at about 0.15 S/m (~1,000 mg/L).

This is why the synthetics exist: to find out when the geophysics can answer a
water-quality question and when it is really answering a question about your
assumptions. A gamma log gives clay CONTENT; spectral gamma or XRD is what
constrains clay TYPE.
""")
