using ArtifactWrappers
using DelimitedFiles
using ClimaLSM
using Plots
using Statistics
# Soil water retention curve data is from
# Wood, Jeffrey D, Gu, Lianhong, Hanson, Paul J, Frankenberg, Christian, & Sack,Lawren. (2022). Supporting biophysical data for "The ecosystem wilting point defines drought response and recovery of a Quercus-Carya forest" [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7477879
af = ArtifactFile(
    url = "https://caltech.box.com/shared/static/8aw3c3zygc7knw94py79826a9ai8hvkb.csv",
    filename = "MOFLUX_SWRC_data.csv",
)
dataset = ArtifactWrapper(@__DIR__, "ozark_swr_data", ArtifactFile[af]);
dataset_path = get_data_folder(dataset);
data = joinpath(dataset_path, "MOFLUX_SWRC_data.csv");
swrc_data = readdlm(data, ','; skipstart = 1)
vwc = swrc_data[4:end, 1] ./ 100 # to convert from percent to a decimal
ψ = swrc_data[4:end, 2] .* 1e6 ./ (1e3 * 9.8) # to convert from MPa to m

# Assume value for this:
θ_r = 0.04
α = 10.0 .^ (log10(5e-3):0.1:1)
n = Array(1.1:0.05:3.0)
porosity = [0.45, 0.5, 0.55]
L2 = zeros(Float64, (length(α), length(n)))

for ν in porosity
    L2 .= 0.0
    for i in 1:length(α)
        for j in 1:length(n)
            hcm = ClimaLSM.Soil.vanGenuchten(; α = α[i], n = n[j])
            S = ClimaLSM.Soil.effective_saturation.(ν, vwc, θ_r)
            ψ_model = ClimaLSM.Soil.matric_potential.(hcm, S)
            L2[i, j] = min(500, sqrt(mean((ψ_model .- ψ) .^ 2)))
        end
    end
    Plots.heatmap(log10.(L2))
    min_idx = argmin(L2)
    Plots.scatter!([min_idx[2]], [min_idx[1]])
    Plots.savefig("./heatmap_ν_$ν.png")
    @info ν, α[min_idx[1]], n[min_idx[2]], minimum(L2)
end
#=
[ Info: (0.45, 0.03154786722400966, 2.1, 44.40710028534568)
[ Info: (0.5, 0.03971641173621406, 2.05, 44.26024841595485)
[ Info: (0.55, 0.049999999999999996, 2.0, 44.34184235590622)
=#
# All very similar, try:
ν = 0.5
α = 0.04
n = 2.05
hcm = ClimaLSM.Soil.vanGenuchten(; α = α, n = n)
θ = 0.1:1e-3:ν
S = ClimaLSM.Soil.effective_saturation.(ν, θ, θ_r)
ψ_model = ClimaLSM.Soil.matric_potential.(hcm, S)
Plots.scatter(-ψ, vwc, label = "Data")
Plots.plot!(-ψ_model, θ, label = "Model")
Plots.plot!(xlabel = "-ψ (m)", ylabel = "θ")
