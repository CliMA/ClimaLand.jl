# Blog post draft: "A sponge, a straw, and fifty millimeters of rain"

Draft of the ClimaLand.jl post for the CliMA software-stack series.

- `blog.md` / `blog.html`: the text (the HTML is the styled preview).
- `snippet.jl`: the reader-facing script printed in the post (bare loam + forest, ~1 min).
- `rain_experiment.jl`: all cases (`CASES`, `NDAYS`, `DT`, `PLANT_A` env vars), writes `results.jls`.
- `plot_results.jl`: figures, water budget table, and the animation from `results.jls`.
- `figs/`, `where_does_the_rain_go.png`: outputs used in the post.

Run from the repo root with `julia +1.12 --startup-file=no --project=.buildkite <script>`.
