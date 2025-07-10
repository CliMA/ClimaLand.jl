using ClimaLand

climaland_dir = pkgdir(ClimaLand)


if length(ARGS) < 1
    error("Must provide site ID on command line")
end

include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/any_run_fluxnet.jl",
    ),
)

# simulated output
run_single_site(ARGS[1], nothing, true)
