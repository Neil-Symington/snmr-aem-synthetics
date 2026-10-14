# SNMR + AEM synthetics

Build a synthetic groundwater section from hydrogeological intervals, derive
what surface NMR and airborne EM would measure over it, and test whether the
inversions can recover the aquifer's architecture, properties and water quality
when the answer is known.

You describe the section as a list of units — each with a porosity, saturation,
clay fraction and pore-water conductivity — and two rock-physics relations give
the two observables:

```
θ_SMR  = φ · S · (1 − clay)^p_clay          free water, what NMR sees
σ_bulk = σ_w · θ_total^m + σ_clay · clay    conductivity, what AEM sees
```

The second relation is the reason this repo exists. One measurement, two
unknowns: a conductive layer may be salty or merely clayey, and EM alone cannot
tell you which. NMR's water content constrains clay from above, because clay
binds water and bound water is invisible to NMR. Whether that is enough to
recover salinity is the question these synthetics are built to answer.

## Quick start

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'     # once, ~5 min

julia --project=. scripts/build_truth.jl config/confined_aquifer.toml
julia --project=. scripts/plot_truth.jl  output/confined_aquifer_truth.png
julia --project=. scripts/forward.jl     output/confined_aquifer_truth.jld2
```

`build_truth.jl` prints the section as built, with the derived observables per
unit measured clear of the gradational boundaries, and writes a JLD2.

## Defining a scenario

A scenario is one TOML file. Units are listed in order of increasing depth and
each runs down to the next one's top; the last extends to `grid.zmax`.

```toml
name = "my_section"

[grid]
dz   = 0.25
zmax = 200.0

[petrophysics]
p_clay     = 1.5    # bound-water exponent
m_archie   = 1.7    # cementation exponent
sigma_clay = 0.5    # S/m per unit clay fraction (surface conduction)

[variability]
seed           = 20261012
corr_length    = 6.0     # m, correlation length of the porosity/clay jitter
log_sd_phi     = 0.10
log_sd_clay    = 0.15
transition_min = 3.0     # gradational boundary thickness, drawn U(min,max)
transition_max = 6.0

[watertable]             # optional; omit for a saturated column
depth      = 6.0
fringe     = 0.8
s_residual = 0.10

[[unit]]
top = 0.0
name = "sand, unsaturated"
phi = 0.33
clay = 0.06
sigma_w = 0.03
# ... more units
```

Two worked examples are in `config/`: a **confined** aquifer enclosed by clay
above and below, and an **unconfined** aquifer that is the top saturated unit.
They differ mainly in where the target sits relative to the conductors, which
turns out to matter a great deal.

### The constants are not innocent

`sigma_clay` is an effective surface conductivity per unit clay fraction, and it
stands in for cation exchange capacity. It varies by close to an order of
magnitude with mineralogy — roughly 0.15 for kaolinite, 0.5 for illite, 2.0 for
smectite. On a fresh aquifer the clay term can be 80–90 % of the total
conductivity, so the salinity you recover is sensitive to a number you probably
do not know. Vary it and see. `p_clay` carries the same dependence and is worth
varying too.

## What is here

| | |
|---|---|
| `src/SNMRAEMSynthetics.jl` | the column builder: config → gradational boundaries → correlated variability → rock physics |
| `src/aem_setup.jl` | SkyTEM-class AEM geometry and AusAEM-standard noise, via HiQGA |
| `scripts/build_truth.jl` | TOML → truth column → JLD2 |
| `scripts/forward.jl` | AEM and NMR forward responses with realistic noise |
| `scripts/plot_truth.jl` | lithology, clay, free water and conductivity on a common depth axis |
| `config/` | two worked scenarios |

## Inversions

Everything above runs on public dependencies. The **AEM** inversion is HiQGA's
and needs nothing further.

The **surface NMR** inversion — trans-dimensional, with the hierarchical noise
model that learns the noise level rather than assuming it — lives in
[`SMRPInversion.jl`](https://github.com/Neil-Symington/SMRPInversion.jl):

```julia
using Pkg
Pkg.add(url = "https://github.com/Neil-Symington/SMRPInversion.jl",
        rev = "add-occam-deterministic-inversion")
```

> **Not yet usable.** As of writing, `src/SMRPInversion.jl` on that branch does
> `include("HierarchicalNoise.jl")` for a file that has never been committed, so
> a fresh clone fails to load. `HierarchicalNoise.jl` and `BootstrapInversion.jl`
> are both untracked locally. They need to be committed and pushed before this
> step works for anyone else.

## A known trap

`makenoisydata!` in HiQGA takes its **own** `rseed` keyword, defaulting to 123,
and reseeds internally — so a `Random.seed!` before calling it does nothing, and
every "different" seed produces identical data. `src/aem_setup.jl` passes `rseed`
through explicitly. Worth knowing if you adapt this code.
