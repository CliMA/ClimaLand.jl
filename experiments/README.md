## ClimaLSM Experiments

These experiment scripts are primarily intended for testing models via
simulation, often comparing to other model output or to data. These run
as part of CI on the Caltech cluster, and must run without error in order
to pass. Many create plots and these plots should be visually assessed before
merging a PR. These experiments do not take the place of unit tests or
tutorials.


We recommend you to use Julia1.9, because the Manifest.toml which is
checked in was created with Julia1.9. If you require using Julia1.8, you
will need to clone the repo, and then run, from ClimaLSM/experiments/:


```
rm Manifest.toml
julia --project
]instantiate
]dev ..
```