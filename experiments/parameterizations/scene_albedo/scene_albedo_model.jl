using DelimitedFiles
using Flux
using MLUtils
using ProgressMeter
using Statistics
using CairoMakie
outpath = "experiments/parameterizations/scene_albedo"
include("experiments/parameterizations/models.jl")
FT = Float32
(net_features, feature_header) = readdlm("experiments/parameterizations/net_features.csv", ',', header=true);
data = net_features;

target_id = [23]
feature_ids = [1,6,7,8,9,10,11,12,16, 17,18,19,20,21,22]
nfeatures = length(feature_ids)
model = Chain(
    Dense(nfeatures => nfeatures, elu),
    Dense(nfeatures => nfeatures, elu),
    Dense(nfeatures => 1, sigmoid),
)
(train, test) = splitobs(transpose(data), at=0.7, shuffle= true);
input = Float32.(train[feature_ids,:]);
truth = Float32.(train[target_id,:]);
test_truth = FT.(test[target_id,:]);
test_input = FT.(test[feature_ids,:]);

loader = Flux.DataLoader((input, truth), batchsize=256, shuffle=true);

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

            Flux.mae(y_hat, y)
        end
        Flux.update!(opt_state, model, grads[1])
        push!(losses_train, loss)  # logging, outside gradient context
    end
    y_hat_test = model(test_input)

    loss_test = Flux.mae(y_hat_test, test_truth)
    push!(losses_test, loss_test)  # logging, outside gradient context
end
# check overfitting
fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1,1], title = "Learning Curves", yscale = log10)
lines!(ax,(1:1:length(losses_train)) ./ length(losses_train),  losses_train, label = "Train")
lines!(ax, (1:1:length(losses_test)) ./ length(losses_test),  losses_test, label = "Test")
CairoMakie.save("experiments/parameterizations/scene_albedo/learning_curves.png", fig)

# Visualize predictions (blend of train and test!)
m = 100000
random_indices = Int.(round.(rand(m).* size(data)[1]))
pred = model(transpose(Float32.(data[random_indices,feature_ids])))
truth = transpose(data[random_indices,target_id])
poor = abs.(pred .- truth) .>0.1;
poor_count = mean(poor)
@show sqrt(mean((pred .- truth).^2))

(metadata, meta_header) = readdlm("metadata.csv", ',', header=true);
metadata = transpose(metadata[random_indices,:])
lat = metadata[1, :]
lon = metadata[2,:]
latg, long, pred_grid = put_on_grid(lat, lon, pred[:])
latg, long, truth_grid = put_on_grid(lat, lon, truth[:])

clims = extrema(pred[.! isnan.(pred)])
fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1,1], title = "Learned")
hm = CairoMakie.heatmap!(ax, latg, long, pred_grid, colorrange = clims)
Colorbar(fig[1, 2], hm)
ax = CairoMakie.Axis(fig[2,1], title = "Bias: Learned - ERA5")
Δ = pred_grid .- truth_grid;
clims = (quantile(Δ[.! isnan.(Δ)],0.01), quantile(Δ[.! isnan.(Δ)],0.99))

hm2 = CairoMakie.heatmap!(ax, latg, long,Δ, colorrange = clims)
Colorbar(fig[2, 2], hm2)
CairoMakie.save("experiments/parameterizations/scene_albedo/albedo_predictions.png", fig)
@show "Learned RMSE:", sqrt(mean((pred[:] .- truth[:]).^2))
@show "Learned bias:", mean(pred[:] .- truth[:])
