using DelimitedFiles
using Flux
using MLUtils
using ProgressMeter
using Statistics
using Plots
include("experiments/parameterizations/models.jl")
FT = Float32
(net_features, feature_header) = readdlm("net_features.csv", ',', header=true);
(nonnormed_input, input_header) = readdlm("nonnormed_input.csv", ',', header=true);
data = hcat(nonnormed_input, net_features);

target_id = [21+4]
LAI_id = [2]
θ_id = [1]
frac_diff_id = [3]
μ_id = [4]
soil_ids = [Array((4+2):(4+12))...]
canopy_ids = [Array((4+13):(4+18))..., 4+1]

(train, test) = splitobs(transpose(data), at=0.7, shuffle= true);
input = Float32.(train[1:24,:]);
truth = Float32.(train[target_id,:]);
test_truth = FT.(test[target_id,:]);
test_input = FT.(test[1:24,:]);

loader = Flux.DataLoader((input, truth), batchsize=256, shuffle=true);

nfeatures_soil = length(soil_ids)
nfeatures_canopy = length(canopy_ids)
model = TotalAlbedo(nfeatures_soil, nfeatures_canopy)
model.logΔα.layers[1].weight[1] = FT(-2)
model.logΔα.layers[1].bias[1] = FT(-1.22)
Flux.trainable(model)
opt = OptimiserChain(ClipGrad(1e-3), Flux.Adam(Float32(5e-5)))
opt_state = Flux.setup(opt, model)
losses_train = []
losses_test = []
@showprogress for epoch in 1:20
    for xy_cpu in loader
        # Unpack batch of data, and move to GPU:
        x, y = xy_cpu
        loss, grads = Flux.withgradient(model) do m
            # Evaluate model and loss inside gradient context:
            y_hat = m(x; soil_ids = soil_ids, canopy_ids=canopy_ids, μ_id=μ_id, frac_diff_id=frac_diff_id, LAI_id= LAI_id, θ_id = θ_id)

            Flux.mae(y_hat, y)
        end
        Flux.update!(opt_state, model, grads[1])
        push!(losses_train, loss)  # logging, outside gradient context
    end
    y_hat_test = model(test_input; soil_ids = soil_ids, canopy_ids=canopy_ids, μ_id=μ_id, frac_diff_id=frac_diff_id, LAI_id= LAI_id, θ_id = θ_id)

    loss_test = Flux.mae(y_hat_test, test_truth)
    push!(losses_test, loss_test)  # logging, outside gradient context
end
# check overfitting
Plots.plot((1:1:length(losses_train)) ./ length(losses_train),  losses_train, yaxis = :log10)
Plots.plot!((1:1:length(losses_test)) ./ length(losses_test),  losses_test)


random_indices = Int.(round.(rand(10000).* size(data)[1]))
pred = model(transpose(Float32.(data[random_indices,:])); soil_ids = soil_ids, canopy_ids=canopy_ids, μ_id=μ_id, frac_diff_id=frac_diff_id, LAI_id= LAI_id, θ_id = θ_id)
truth = transpose(data[random_indices,target_id])
poor = abs.(pred .- truth) .>0.1;
poor_count = mean(poor)
@show sqrt(mean((pred .- truth).^2))
(metadata, meta_header) = readdlm("metadata.csv", ',', header=true);
metadata = transpose(metadata[random_indices,:])
lat = metadata[1, :]
lon = metadata[2,:]

dry_albedo = model.dry_0(transpose(Float32.(data[random_indices,:]))[soil_ids, :])
wet_albedo = model.wet_0(transpose(Float32.(data[random_indices,:]))[soil_ids, :])
latg, long, xgrid = put_on_grid(lat, lon, wet_albedo[:] .- dry_albedo[:])



(metadata, meta_header) = readdlm("metadata.csv", ',', header=true);
metadata = transpose(metadata[random_indices,:])
Plots.scatter(metadata[2, poor[:]], metadata[1, poor[:]])

# clm prediction
(clm_features, clm_header) = readdlm("clm_albedos.csv", ',', header=true);
clm_features = transpose(clm_features[random_indices,:])

α_leaf = (clm_features[1,:] .+ clm_features[2,:]) ./ 2
τ_leaf = (clm_features[3,:] .+ clm_features[4,:]) ./ 2
θ = transpose(Float32.(data[random_indices,θ_id]))[:]
Δ = @. max(0.11 − 0.40*θ, 0)
α_soil_PAR = @. min(clm_features[5,:] .+ Δ, clm_features[6,:])
α_soil_NIR = @. min(clm_features[7,:] .+ Δ, clm_features[8,:])
α_soil = (α_soil_PAR .+ α_soil_NIR) ./ 2

LAI = data[random_indices,LAI_id[1]]
μ = data[random_indices,μ_id[1]]
frac_diff = data[random_indices,frac_diff_id[1]]
clm_pred = simpler_canopy_sw_rt_two_stream.(α_leaf, τ_leaf, LAI, μ, α_soil, frac_diff)
@show sqrt(mean((clm_pred[:] .- truth[:]).^2))

slice =1:1:1000
Plots.scatter(pred[slice], pred[:][slice] .- truth[:][slice], label = "new model")
#Plots.scatter!(μ[slice], clm_pred[:][slice] .- truth[:][slice], label = "clm")
