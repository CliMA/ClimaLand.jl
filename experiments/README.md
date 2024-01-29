## ClimaLand Experiments

These experiment scripts are primarily intended for testing models via
simulation, often comparing to other model output or to data. These run
as part of CI on the Caltech cluster, and must run without error in order
to pass. Many create plots and these plots should be visually assessed before
merging a PR. These experiments do not take the place of unit tests or
tutorials.


We recommend you to use Julia1.10, because the Manifest.toml which is
checked in was created with Julia1.10. If you require using an earlier
version of Julia, you will need to clone the repo,
and then run, from ClimaLand/experiments/:


```
rm Manifest.toml
julia --project
]instantiate
]dev ..
```
