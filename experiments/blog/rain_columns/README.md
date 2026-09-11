# Blog post draft: "A sponge, a straw, and fifty millimeters of rain"

Draft of the ClimaLand.jl post for the CliMA software-stack series.

- `blog.md`: the post. `blog.html` is the styled preview, rendered from it by `build_html.py`
  (which also inserts the code excerpt from `snippet.jl`); `style.css` is the preview's stylesheet.
- `snippet.jl`: the reader-facing script linked from the post (bare loam + forest, ~1 min).
- `rain_experiment.jl`: all five cases; env vars `CASES`, `NDAYS`, `DT`, `PLANT_A`, `STORM_MM`, `OUTDIR`.
- `plot_results.jl`: figures, water budget table, and the animation.
- `figs/`, `where_does_the_rain_go.png`: the outputs used in the post (committed).

To regenerate everything, from the repository root:

```
julia +1.12 --startup-file=no --project=.buildkite experiments/blog/rain_columns/rain_experiment.jl   # OUTDIR=out (default)
STORM_MM=0 OUTDIR=experiments/blog/rain_columns/out_control julia +1.12 --startup-file=no --project=.buildkite experiments/blog/rain_columns/rain_experiment.jl
cd experiments/blog/rain_columns
julia +1.12 --startup-file=no --project=../../../.buildkite plot_results.jl   # reads out/ and out_control/
julia +1.12 --startup-file=no --project=../../../.buildkite snippet.jl        # writes where_does_the_rain_go.png
python3 build_html.py
```

The budget figure and the numbers in the text are storm run minus no-storm control, so both runs are needed.
