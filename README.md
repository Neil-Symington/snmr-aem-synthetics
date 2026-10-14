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

## Examples first

Three short scripts in `examples/`, none longer than a few seconds and none
needing an inversion:

```bash
julia --project=. examples/01_hello_column.jl      # smallest working section
julia --project=. examples/02_build_in_code.jl     # same, without TOML
julia --project=. examples/03_why_clay_matters.jl  # why conductivity != salinity
```

The third is the one to read. It infers salinity from a *perfect* conductivity
measurement while varying only the assumed clay mineralogy, and the answer moves
by a factor of five. See `examples/README.md`.

## Defining a scenario

A scenario is one TOML file. Units are contiguous and listed top down. Give each
one a `thickness` and the tops accumulate automatically, or give `top` on every
unit, or mix the two — they are checked against each other.

```toml
[[unit]]
name = "target aquifer"
lithology = "sand"        # fills porosity and clay from the preset table
thickness = 22.0          # m
phi = 0.31                # override the preset wherever you know better
sigma_w_mgL = 520         # salinity in mg/L, or sigma_w in S/m
sd_phi = 0.05             # this unit is more uniform than the section default
```

`config/template.toml` shows every option with commentary. Three things worth
knowing:

**Lithology presets** fill porosity and clay so you can write the section in the
terms you actually know. Available: `gravel`, `clean sand`, `sand`, `silty
sand`, `sandy clay`, `clay`, `weathered bedrock`, `saprock`, `fresh bedrock`.
They are indicative starting points, not measurements — anything you state
explicitly wins.

**The truth is jiggled, not piecewise constant.** Each property gets a smooth,
depth-correlated field scaled by a log-normal spread, so a unit varies around
the value you gave rather than being a flat slab. Spreads can be set section
wide in `[variability]`, per lithology (a clean sand is more uniform than a
weathered profile), or per unit. Porosity and clay vary by default; pore-water
conductivity does **not**, so salinity stays exactly as stated unless you set
`log_sd_sigma_w`.

**Boundaries are gradational**, over a thickness drawn uniformly between
`transition_min` and `transition_max`, so no contact is artificially sharp.

### Is your section a fair test?

`build_truth.jl` runs `check_column` and prints a conductivity budget — how much
of each unit's conductivity is carried by pore water and how much by clay
surfaces — then flags the mistakes that look fine on paper:

- a target where clay carries most of the conductivity, with the error
  amplification quoted, so you know before inverting how sensitive the salinity
  answer will be;
- adjacent units that both clear the extractable bar and so will be detected as
  **one** aquifer rather than two;
- a large salinity step across permeable ground with no confining layer;
- pore water freshening with depth;
- a target that is a resistor enclosed by conductors — the case where airborne
  EM resolves the unit poorly and the salinity may be set by the regularisation
  rather than the data.

On `config/confined_aquifer.toml` it independently identifies the two findings
that took a week to establish empirically.

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
| `src/SNMRAEMSynthetics.jl` | the column builder: config → geometry → correlated variability → rock physics, plus `check_column` |
| `src/Lithology.jl` | porosity and clay presets by material, and mg/L ↔ S/m |
| `src/aem_setup.jl` | SkyTEM-class AEM geometry and AusAEM-standard noise, via HiQGA |
| `scripts/build_truth.jl` | TOML → truth column → JLD2 |
| `scripts/forward.jl` | AEM and NMR forward responses with realistic noise |
| `scripts/plot_truth.jl` | lithology, clay, free water and conductivity on a common depth axis |
| `config/` | a minimal scenario, a fully commented template, and the two realistic ones |
| `examples/` | three short runnable demonstrations |

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

Verified to install and load from a clean environment. `HierarchicalNoise.jl`
and `BootstrapInversion.jl` were untracked upstream until October 2026 — the
module included files that had never been committed, so a fresh clone failed to
load. That is fixed on the branch above.

What the hierarchical noise model does: the per-gate variance shape is held
fixed and only its scalar scale is sampled, jointly with the water content, as a
nuisance parameter. The inversion therefore learns its own noise level rather
than being handed one. That matters because an assumed noise level that is too
low manufactures structure and one that is too high erases it — and on field
data you do not know which you have.

## A known trap

`makenoisydata!` in HiQGA takes its **own** `rseed` keyword, defaulting to 123,
and reseeds internally — so a `Random.seed!` before calling it does nothing, and
every "different" seed produces identical data. `src/aem_setup.jl` passes `rseed`
through explicitly. Worth knowing if you adapt this code.
