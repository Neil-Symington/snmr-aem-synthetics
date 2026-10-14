#!/usr/bin/env julia
# Example 1 — the smallest thing that does something.
#
#   julia --project=. examples/01_hello_column.jl
#
# Builds a three-unit section from config/minimal.toml and prints what the two
# instruments would see. Runs in a couple of seconds; no inversion involved.

include(joinpath(@__DIR__, "..", "src", "SNMRAEMSynthetics.jl"))
using .SNMRAEMSynthetics, Printf

units, cfg = load_config(joinpath(@__DIR__, "..", "config", "minimal.toml"))
col = build_column(units, cfg)

println("Three units in, two observables out:\n")
show(stdout, col)

println("""

Reading the two right-hand columns:

  theta_SMR  free water — what surface NMR measures, more or less directly.
             The aquifer is high and the clay is low, because clay binds water
             and bound water is invisible to NMR.

  sigma_b    bulk conductivity — what airborne EM measures. Note that the clay
             aquitard is conductive mostly BECAUSE it is clay, while the
             basement is resistive despite having the saltiest pore water in
             the section, because it has almost no porosity to conduct through.

That second point is the whole difficulty: conductivity mixes clay and salt,
and one measurement cannot separate two unknowns. Example 3 puts numbers on it.
""")
