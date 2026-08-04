# Prognostic Carbon Model

The prognostic carbon model carries live vegetation carbon as four prognostic
pools and predicts biomass from the carbon the canopy fixes, rather than
prescribing it. It is a *wrapper* around an existing LAI model: it occupies the
same `:biomass` component slot, forwards everything the wrapped model provides,
and adds the pools on top.

## Pools and fluxes

Four pools, all in kg C m⁻² of ground:

| pool | contents |
|---|---|
| ``C_{sugar}`` | non-structural carbon; the buffer GPP enters and allocation leaves |
| ``C_{leaf}`` | leaf carbon |
| ``C_{stem}`` | woody carbon, including sapwood |
| ``C_{root}`` | fine-root carbon |

Gross primary production enters the sugar pool, maintenance respiration and
allocation leave it, and the three structural pools each gain a fraction of what
is allocated and lose carbon to litter at their own turnover time:

```math
\begin{align}
\frac{dC_{sugar}}{dt} &= \text{GPP} - R_m - S \\
\frac{dC_{leaf}}{dt}  &= a f_{leaf} S - C_{leaf}/\tau_{leaf} \\
\frac{dC_{stem}}{dt}  &= a f_{stem} S - C_{stem}/\tau_{stem} \\
\frac{dC_{root}}{dt}  &= a f_{root} S - C_{root}/\tau_{root}
\end{align}
```

``S`` is the allocation flux out of the sugar pool, ``a`` is the construction
efficiency, so ``(1-a)S`` is growth respiration, and the ``f`` are the allocation
fractions described below. Litter leaving the three structural pools is passed to
the soil carbon pool by the integrated model, closing the balance against the
microbial drawdown.

Allocation is throttled by how full the sugar pool is relative to a target set by
the standing structural biomass, so a canopy cannot allocate carbon it has not
fixed, and maintenance respiration shuts down as the sugar pool empties rather
than driving it negative.

## Allocation fractions are global constants, modified by climate

The three fractions sum to one at every column and every timestep. There are no
plant functional types anywhere in the model: a PFT map is a boundary condition
that cannot be extrapolated to a climate that has not been observed, which
defeats the purpose of a land model inside a climate model. Where allocation must
vary in space it varies with **climate**, through two dimensionless multipliers on
the woody fraction.

### Precipitation: the woody fraction

Mean annual precipitation is the classical first-order control on maximum woody
cover in seasonally dry systems (Sankaran et al., 2005). Woody allocation follows
a saturating ramp in it:

```math
w(P) = \frac{x^n}{1 + x^n}, \qquad x = P_{annual}/P_{half}
```

so a dry column builds almost no wood and a wet one approaches the full
``f_{stem}``.

### Rainfall seasonality: the water-deficit limit

Mean annual precipitation is an annual *total*, and it cannot separate wet
savanna from wet forest — the Cerrado, the Sahel fringe, the Pampas and miombo
all receive forest-sized rainfall and carry no forest. What differs is how that
rain is distributed through the year: 1500 mm spread evenly supports closed
canopy, while the same 1500 mm delivered in five months, followed by a seven-month
drought that dries fuel and limits tree establishment, supports savanna.

The model therefore accumulates an **annual water deficit** ``D`` — the shortfall
of monthly precipitation below evaporative demand, summed over the year — and
applies a second multiplier:

```math
\ell(D) = \frac{1}{1 + (D/D_{half})^{n_D}}
```

Evaporative demand carries temperature, as a ramp from a floor at freezing to a
full reference at 20 °C above it:

```math
E(T) = E_{ref}\left[\phi + (1-\phi)\,\mathrm{clamp}\!\left(\tfrac{T - 273.15}{20},\,0,\,1\right)\right]
```

The floor ``\phi`` is not cosmetic. A *constant* reference assumes tropical
demand everywhere and scores a cold, dry boreal cell as severely water-stressed
when the reason it is dry is that it is frozen, which suppresses boreal forest
where observations show closed canopy. A pure ramp to zero protects the boreal
but gives up most of the signal. The floor interpolates between the two.

Stem allocation is the product of both limits, and whatever stem gives up goes to
root, so the fractions still sum to one:

```math
f_{stem} = f_{stem,0}\; w(P_{annual})\; \ell(D_{annual}), \qquad
f_{root} = 1 - f_{leaf} - f_{stem}
```

!!! note "The averaging window is load-bearing"
    The deficit is taken against a 30-day running precipitation total, not the
    instantaneous rate. Precipitation is zero at almost every instant, so a
    deficit taken against the instantaneous rate collapses to the annual
    evaporative demand in every column and carries no seasonality information at
    all — the mechanism would appear implemented and measure nothing.

### Temperature: stem turnover

Stem longevity lengthens in the cold, so ``\tau_{stem}`` is scaled by a
``Q_{10}``-style factor in mean annual temperature, capped so that a very cold
column cannot acquire an unbounded turnover time.

## Time-integrated climate state

The climate quantities above are prognostic, not read from a map. The model
declares its own running means and sums, so it works under a prescribed-LAI
configuration where the optimal-LAI model and its own integrals are absent:

| variable | reduction | meaning |
|---|---|---|
| ``T_{annual}`` | running mean | mean annual temperature |
| ``P_{annual}`` | 1-year running sum | mean annual precipitation |
| ``P_{month}`` | 30-day running sum | monthly precipitation total |
| ``T_{month}`` | 30-day running mean | monthly mean temperature |
| ``D_{annual}`` | 1-year running sum | annual water deficit |

!!! warning "Spin-up"
    ``D_{annual}`` has a memory of order two years and feeds a pool whose
    turnover is decades. A coupled run shorter than several times that memory
    **under-suppresses woody allocation by construction** — the deficit has not
    finished accumulating. Initial conditions assume an aseasonal climate, so the
    error is in the direction of too little suppression, which shows up against
    observations rather than hiding.

## Coupling to GPP and LAI

The coupling is currently **one-way**: the pools consume GPP and LAI and feed
nothing back. This is a deliberate staging choice, and it makes an important
property structural rather than a matter of care — adding the carbon model cannot
change GPP or LAI at all, so a simulation's photosynthesis and canopy structure
are bit-identical with the pools on and off. That is verified directly by the
`Carbon model does not change LAI or GPP` testset and, at scale, by a
twenty-site battery.

The one sanctioned exception is autotrophic respiration:
[`PoolBasedAutotrophicRespirationModel`](@ref) can take ``R_a = R_m + R_g`` from
the pools instead of the default parameterization, since a model carrying
explicit sapwood and root pools can compute maintenance respiration from them.

## Known limitations

These are properties of the current formulation, not bugs:

- **Disturbance is absent.** The model builds forest wherever climate permits it.
  The residual overprediction is concentrated in a small minority of treeless
  cells whose rainfall *is* forest rainfall; the water-deficit limit reaches many
  but not all of them. Fire, grazing and cultivation are not represented, and a
  prescribed burned-area map was deliberately rejected as being as
  unextrapolatable as a PFT map.
- **Root allocation is not observationally constrained.** No gridded product
  constrained the belowground fractions, so they are global constants that the
  biomass comparison cannot test.
- **The deficit depends on the forcing's rainfall seasonality**, which is not
  uniformly good. Where a forcing dataset smooths the dry season, the limit
  weakens accordingly.

## API

```@docs
ClimaLand.Canopy.PrognosticCarbonModel
ClimaLand.Canopy.PrognosticCarbonParameters
ClimaLand.Canopy.woody_fraction
ClimaLand.Canopy.seasonality_limit
ClimaLand.Canopy.monthly_pet
ClimaLand.Canopy.tau_stem_scale
```
