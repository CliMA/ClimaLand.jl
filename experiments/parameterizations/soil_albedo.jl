# Try making heatmaps of variables - are they reasonable? do they correlate with soil type?
# Look at CLM prediction - is moisture variable really wrong?

import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaUtilities
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
using ClimaCore
using Plots

import ClimaParams as CP
using Dates
using Flux
using MLUtils
using ProgressMeter
using Statistics
using Plots

using ClimaLand
import ClimaLand.Parameters as LP
start_date = DateTime(2008)
time_interpolation_method =LinearInterpolation()
regridder_type = :InterpolationsRegridder
nelements = (101, 15)
FT = Float32
domain = ClimaLand.global_domain(FT; nelements)
surface_space = domain.space.surface
subsurface_space = domain.space.subsurface

spatially_varying_soil_params =
    ClimaLand.default_spatially_varying_soil_parameters(
        subsurface_space,
        surface_space,
        FT,
    )

(;
 ν,
 ν_ss_om,
 ν_ss_quartz,
 ν_ss_gravel,
 hydrology_cm,
 K_sat,
 θ_r,
 ) = spatially_varying_soil_params

modis_lai_artifact_path =
    ClimaLand.Artifacts.modis_lai_forcing_data_path()
modis_lai_ncdata_path =
    joinpath(modis_lai_artifact_path, "Yuan_et_al_2008_1x1.nc")

LAIfunction = ClimaLand.prescribed_lai_modis(
    modis_lai_ncdata_path,
    surface_space,
    start_date;
    time_interpolation_method = time_interpolation_method,
)


era5_path = "../ClimaArtifacts/era5_soil_albedo_data/soil_a_data_2008_1.0x1.0_avg.nc";
function α(net, down)
    if down < 0.1
        return 0.0
    else
        return 1 - net/down
    end
end

albedo = TimeVaryingInput(
    era5_path,
    ["avg_snswrf","avg_sdswrf"],
    surface_space;
    start_date,
    regridder_type,
    method = time_interpolation_method,
    compose_function = (net, down) -> α.(net, down)
)
era5_path = "../ClimaArtifacts/era5_soil_albedo_data/soil_a_data_2008_1.0x1.0_inst.nc";
snow_present = TimeVaryingInput(
    era5_path,
    "sd",
    surface_space;
    start_date,
    regridder_type,
    method = time_interpolation_method,
    compose_function = (sd) -> sd .> 0.001
)

θ = TimeVaryingInput(
    era5_path,
    "swvl1",
    surface_space;
    start_date,
    regridder_type,
    method = time_interpolation_method,
)


lai_field = ClimaCore.Fields.zeros(surface_space)
albedo_field = ClimaCore.Fields.zeros(surface_space)
θ_field = ClimaCore.Fields.zeros(surface_space)
snow_field = ClimaCore.Fields.zeros(surface_space)
mask = ClimaCore.Fields.zeros(surface_space)
N = 0
for (i,t) in enumerate(θ.data_handler.available_dates)
    evaluate!(lai_field, LAIfunction, t)
    evaluate!(albedo_field, albedo, t)
    evaluate!(θ_field, θ, t)
    evaluate!(snow_field, snow_present, t)
    @. mask = ! ((snow_field ≈ 1.0) || (lai_field > 0.001) || (albedo_field < 0.01))
    array_mask = parent(mask)[:] .≈ 1
    N_i = sum(array_mask)
    @show i
    @show N_i
    N += sum(array_mask)
end
data = zeros(N, 9)
id = 1
for (i,t) in enumerate(θ.data_handler.available_dates)
    evaluate!(lai_field, LAIfunction, t)
    evaluate!(albedo_field, albedo, t)
    evaluate!(θ_field, θ, t)
    evaluate!(snow_field, snow_present, t)
    @. mask = ! ((snow_field ≈ 1.0) || (lai_field > 0.001) || (albedo_field < 0.01))
    array_mask = parent(mask)[:] .≈ 1
    albedo_i = parent(albedo_field)[:][array_mask]
    @assert (minimum(albedo_i)) >= 0.01
    ν_i = parent(ν)[end, :,:,:,:][array_mask]
    ν_ss_om_i = parent(ν_ss_om)[end, :,:,:,:][array_mask]
    ν_ss_quartz_i = parent(ν_ss_quartz)[end,:,:,:,:][array_mask]
    ν_ss_gravel_i = parent(ν_ss_gravel)[end,:,:,:,:][array_mask]
    vg_α_i = parent(hydrology_cm.α)[end,:,:,:,:][array_mask]
    vg_n_i = parent(hydrology_cm.n)[end,:,:,:,:][array_mask]
    K_sat_i = parent(K_sat)[end,:,:,:,:][array_mask]
    θ_i = parent(θ_field)[:][array_mask]
    albedo_i = parent(albedo_field)[:][array_mask]
    N_i =  sum(array_mask)
    data[id:id+N_i-1, 1] .= θ_i # Why is this so small?
    data[id:id+N_i-1, 2] .= ν_ss_om_i
    data[id:id+N_i-1, 3] .= ν_ss_quartz_i
    data[id:id+N_i-1, 4] .= ν_ss_gravel_i
    data[id:id+N_i-1, 5] .= 1.0 ./ vg_α_i
    data[id:id+N_i-1, 6] .= vg_n_i .- 1
    data[id:id+N_i-1, 7] .= log10.(K_sat_i)
    data[id:id+N_i-1, 8] .= ν_i
    data[id:id+N_i-1, 9] .= albedo_i
    id = id+N_i
end



# Train a model
# zmuv
data_normed = copy(data);
nfeatures = 8
target_id = [nfeatures+1]
data_normed[:,1:nfeatures] .= (data[:, 1:nfeatures] .- mean(data[:, 1:nfeatures], dims = 1)) ./ std(data[:, 1:nfeatures], dims = 1)
# train test split]
ndata = length(data[:,1])
(train, test) = splitobs(transpose(data_normed), at=0.7, shuffle= true);
input = Float32.(train[1:nfeatures,:])
truth = Float32.(train[target_id,:])
test_truth = FT.(test[target_id,:])
test_input = FT.(test[1:nfeatures,:])
loader = Flux.DataLoader((input, truth), batchsize=128, shuffle=true);


model = Chain(
    Dense(8 => 16, elu),      # activation function inside layer
    BatchNorm(16),
    Dense(16 => 16, elu),      # activation function inside layer
    BatchNorm(16),
    Dense(16 => 1, sigmoid))


opt = OptimiserChain(ClipGrad(1e-3), Flux.Adam(Float32(5e-5)))
opt_state = Flux.setup(opt, model)
losses_train = []
losses_test = []
@showprogress for epoch in 1:30
    for xy_cpu in loader
        # Unpack batch of data, and move to GPU:
        x, y = xy_cpu
        loss, grads = Flux.withgradient(model) do m
            # Evaluate model and loss inside gradient context:
            y_hat = m(x)
            Flux.mse(y_hat, y)
        end
        Flux.update!(opt_state, model, grads[1])
        push!(losses_train, loss)  # logging, outside gradient context
    end
    loss_test = Flux.mse(model(test_input), test_truth)
    push!(losses_test, loss_test)  # logging, outside gradient context
end
# check overfitting
Plots.plot((1:1:length(losses_train)) ./ length(losses_train),  losses_train, yaxis = :log10)
Plots.plot!((1:1:length(losses_test)) ./ length(losses_test),  losses_test)

train_pred = model(input)
poor_count_train = mean(abs.(train_pred .- truth) .>0.05)
Plots.scatter(truth[:], abs.(train_pred[:] .- truth[:]))
