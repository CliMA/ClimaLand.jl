using DelimitedFiles
using Flux
using MLUtils
using ProgressMeter
using Statistics
using CairoMakie
using CurveFit
outpath = "experiments/parameterizations/no_tau"
include("experiments/parameterizations/models.jl")
FT = Float32
(net_features, feature_header) = readdlm("experiments/parameterizations/net_features.csv", ',', header=true);
(nonnormed_input, input_header) = readdlm("experiments/parameterizations/nonnormed_input.csv", ',', header=true);
(metadata, meta_header) = readdlm("experiments/parameterizations/metadata.csv", ',', header=true);
(clm_features, clm_header) = readdlm("experiments/parameterizations/clm_albedos.csv", ',', header=true);



data = hcat(nonnormed_input, net_features);
# train test split]
target_id = [23+4]
LAI_id = [2]
θ_id = [1]
frac_diff_id = [3]
μ_id = [4]
soil_ids = [Array((4+6):(4+11))...,4+1]
canopy_ids = [Array((4+17):(4+21))..., 4+1]
feature_ids = [soil_ids..., canopy_ids...]
(train, test) = splitobs(transpose(data), at=0.7, shuffle= true);
input = Float32.(train[1:26,:]);
truth = Float32.(train[target_id,:]);
test_truth = FT.(test[target_id,:]);
test_input = FT.(test[1:26,:]);

loader = Flux.DataLoader((input, truth), batchsize=256, shuffle=true);

nfeatures_soil = length(soil_ids)
nfeatures_canopy = length(canopy_ids)
model = TotalAlbedoSingleSoilNoTau(nfeatures_soil, nfeatures_canopy)
model.logΔα.layers[1].weight[1] = FT(-2)
model.logΔα.layers[1].bias[1] = FT(-1.22)
model.LAI_sigmoid.layers[1].weight[1] = FT(3.0)
model.LAI_sigmoid.layers[1].bias[1] = FT(-1.0)
Flux.trainable(model)
opt = OptimiserChain(ClipGrad(1e-3), Flux.Adam(Float32(8e-5)))
opt_state = Flux.setup(opt, model)
losses_train = []
losses_test = []
@showprogress for epoch in 1:30
    for xy_cpu in loader
        # Unpack batch of data, and move to GPU:
        x, y = xy_cpu
        loss, grads = Flux.withgradient(model) do m
            # Evaluate model and loss inside gradient context:
            y_hat = m(x; soil_ids = soil_ids, canopy_ids=canopy_ids, μ_id=μ_id, frac_diff_id=frac_diff_id, LAI_id= LAI_id)

            Flux.mae(y_hat, y)
        end
        Flux.update!(opt_state, model, grads[1])
        push!(losses_train, loss)  # logging, outside gradient context
    end
    y_hat_test = model(test_input; soil_ids = soil_ids, canopy_ids=canopy_ids, μ_id=μ_id, frac_diff_id=frac_diff_id, LAI_id= LAI_id)

    loss_test = Flux.mae(y_hat_test, test_truth)
    push!(losses_test, loss_test)  # logging, outside gradient context
end
# check overfitting
outpath = "experiments/parameterizations/no_tau"
fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1,1], title = "Learning Curves", yscale = log10)
lines!(ax,(1:1:length(losses_train)) ./ length(losses_train),  losses_train, label = "Train")
lines!(ax, (1:1:length(losses_test)) ./ length(losses_test),  losses_test, label = "Test")
CairoMakie.save("$outpath/learning_curves.png", fig)

# Visualize predictions (blend of train and test!)
m = 10000
random_indices = Int.(round.(rand(m).* size(data)[1]))
pred = model(transpose(Float32.(data[random_indices,:])); soil_ids = soil_ids, canopy_ids=canopy_ids, μ_id=μ_id, frac_diff_id=frac_diff_id, LAI_id= LAI_id)
truth = transpose(data[random_indices,target_id])
poor = abs.(pred .- truth) .>0.1;
poor_count = mean(poor)
@show sqrt(mean((pred .- truth).^2))


# clm prediction

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


metadata = transpose(metadata[random_indices,:])
lat = metadata[1, :]
lon = metadata[2,:]



# Soil albedo plots
soil_albedo = model.α_soil(transpose(Float32.(data[random_indices,:]))[soil_ids, :])
latg, long, sa_grid = put_on_grid(lat, lon, soil_albedo[:])
_,_,clm_sa_grid = put_on_grid(lat, lon, α_soil[:])

ratio_α_τ_leaf = model.ratio_α_τ_leaf(transpose(Float32.(data[random_indices,:]))[canopy_ids, :])[:]
sum_α_τ_leaf = model.sum_α_τ_leaf(transpose(Float32.(data[random_indices,:]))[canopy_ids, :])[:]
leaf_tau  = sum_α_τ_leaf ./ (1 .+  ratio_α_τ_leaf)
leaf_albedo = sum_α_τ_leaf .- leaf_tau
latg, long, la_grid = put_on_grid(lat, lon, leaf_albedo[:].* (LAI[:] .> 0.05))
latg, long, lt_grid = put_on_grid(lat, lon, leaf_tau[:] .* (LAI[:] .> 0.05) )
latg, long, clm_la_grid = put_on_grid(lat, lon, α_leaf[:])
latg, long, clm_lt_grid = put_on_grid(lat, lon, τ_leaf[:])
latg, long, pred_grid = put_on_grid(lat, lon, pred[:])
latg, long, clm_pred_grid = put_on_grid(lat, lon, clm_pred[:])
latg, long, truth_grid = put_on_grid(lat, lon, truth[:])


fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1,1], title = "Learned")
hm = CairoMakie.heatmap!(ax, latg, long, sa_grid, colorrange = [0,0.5])
Colorbar(fig[1, 2], hm)
ax = CairoMakie.Axis(fig[2,1], title = "CLM")

hm2 = CairoMakie.heatmap!(ax, latg, long, clm_sa_grid, colorrange = [0,0.5])
Colorbar(fig[2, 2], hm2)
CairoMakie.save("$outpath/soil_albedo.png", fig)


fig = CairoMakie.Figure()
xb, yb, counts = put_on_grid_and_count(θ[:], soil_albedo[:])
(b,m) = CurveFit.linear_fit(θ[θ .> 0.15],soil_albedo[:][θ .> 0.15])
ax = CairoMakie.Axis(fig[1,1], xlabel= "θ/ν", ylabel = "soil albedo", title = "Learned; slope = $m")
hm = CairoMakie.contourf!(ax, xb, yb, counts)
lines!(ax, xb, xb .* m .+b, linestyle = :dash, color = :red)
Colorbar(fig[1, 2], hm)
xb, yb, counts = put_on_grid_and_count(θ[:], α_soil[:])
(b,m) = CurveFit.linear_fit(θ[θ .> 0.15],α_soil[:][θ .> 0.15])
ax = CairoMakie.Axis(fig[2,1], xlabel= "θ/ν", ylabel = "soil albedo",  title = "CLM; slope = $m")
hm2 = CairoMakie.contourf!(ax, xb, yb, counts)
lines!(ax, xb, xb .* m .+b, linestyle = :dash, color = :red)

Colorbar(fig[2, 2], hm)
CairoMakie.save("$outpath/soil_albedo_moisture.png", fig)

fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1,1], title = "Learned")
clims = extrema(la_grid[.! isnan.(la_grid)])
hm = CairoMakie.heatmap!(ax, latg, long, la_grid, colorrange = clims)
Colorbar(fig[1, 2], hm)
ax = CairoMakie.Axis(fig[2,1], title = "CLM")
clims = extrema(clm_la_grid[.! isnan.(clm_la_grid)])
hm2 = CairoMakie.heatmap!(ax, latg, long, clm_la_grid, colorrange = clims)
Colorbar(fig[2, 2], hm2)
CairoMakie.save("$outpath/leaf_albedo.png", fig)


fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1,1], title = "Learned")
clims = extrema(lt_grid[.! isnan.(lt_grid)])
hm = CairoMakie.heatmap!(ax, latg, long, lt_grid, colorrange = clims)
Colorbar(fig[1, 2], hm)
ax = CairoMakie.Axis(fig[2,1], title = "CLM")
clims = extrema(clm_lt_grid[.! isnan.(clm_lt_grid)])
hm2 = CairoMakie.heatmap!(ax, latg, long, clm_lt_grid, colorrange = clims)
Colorbar(fig[2, 2], hm2)
CairoMakie.save("$outpath/leaf_transmissivity.png", fig)


clims = extrema(pred[.! isnan.(pred)])
fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1,1], title = "Learned")
hm = CairoMakie.heatmap!(ax, latg, long, pred_grid, colorrange = clims)
Colorbar(fig[1, 2], hm)
ax = CairoMakie.Axis(fig[2,1], title = "CLM")

hm2 = CairoMakie.heatmap!(ax, latg, long, clm_pred_grid, colorrange = clims)
Colorbar(fig[2, 2], hm2)
CairoMakie.save("$outpath/surface_albedo.png", fig)

fig = CairoMakie.Figure()
latg, long, LAI_grid = put_on_grid(lat, lon, LAI[:])

ax = CairoMakie.Axis(fig[1,1], title = "MODIS")
clims = extrema(LAI_grid[.! isnan.(LAI_grid)])
hm = CairoMakie.heatmap!(ax, latg, long, LAI_grid, colorrange = clims)
Colorbar(fig[1, 2], hm)
CairoMakie.save("$outpath/LAI.png", fig)



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
CairoMakie.save("experiments/parameterizations/two_stream_models/single_soil/albedo_predictions.png", fig)

@show "CLM RMSE:", sqrt(mean((clm_pred[:] .- truth[:]).^2))
@show "Learned RMSE:", sqrt(mean((pred[:] .- truth[:]).^2))
@show "CLM bias:", mean(clm_pred[:] .- truth[:])
@show "Learned bias:", mean(pred[:] .- truth[:])
