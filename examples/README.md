# Examples

Three short scripts, in order. None needs an inversion and each runs in a few
seconds, so they are safe to work through before committing to anything long.

```bash
julia --project=. examples/01_hello_column.jl
julia --project=. examples/02_build_in_code.jl
julia --project=. examples/03_why_clay_matters.jl
```

**`01_hello_column.jl` — the smallest thing that does something.** Builds a
three-unit section (sand / clay / basement) from `config/minimal.toml` and
prints what each instrument would see. Start here to check the install works
and to see the shape of the output.

**`02_build_in_code.jl` — the same without TOML.** The config is just a `Dict`
and units are a plain struct, so sections can be generated in a loop. Ends with
a small sweep over aquifer thickness, which is the pattern you would extend for
a scenario suite or a Monte Carlo.

**`03_why_clay_matters.jl` — the one worth reading.** Splits the aquifer's
conductivity into the part carried by pore water and the part carried by clay
surfaces, then infers salinity from a *perfect* conductivity measurement while
varying only the assumed clay mineralogy. The inferred salinity moves by a
factor of five across the realistic range, and two assumptions are inadmissible
outright. It is the clearest short demonstration of why a conductivity
measurement alone does not give you water quality.

## Then what

```bash
julia --project=. scripts/plot_truth.jl output/minimal_truth.png
julia --project=. scripts/forward.jl    output/minimal_truth.jld2
```

`config/confined_aquifer.toml` and `config/unconfined_aquifer.toml` are the two
realistic sections used in the Murray Basin study — a confined aquifer enclosed
by clay, and an unconfined aquifer as the top saturated unit. They are a good
starting point for a section of your own.
