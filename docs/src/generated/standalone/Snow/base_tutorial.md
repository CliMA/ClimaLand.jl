```@meta
EditURL = "https://github.com/CliMA/ClimaLand.jl/../../../../.."
```
# Seasonal Snow Timeseries Generation with a Neural Network

This tutorial explains how to make use of the code developed
for forecasting seasonal snow depth evolution, using a neural network with
structurally-enforced constraints. The following steps through a
basic use-case of the system on an already-cleaned dataset, though
exploration of optional keyword arguments in the developed code and
additional tools for scraping data (explained in the [data tutorial](../data_tutorial/))
provide for a richer set of functionality.

The updates of the neural snow model follow the equation

``
\frac{dz}{dt} = \mathcal{M}\left(z, SWE, φ, R, v, T_{air}, P_{snow}\right),
``

where

- ``t`` is the time (s),

- ``z`` is the snow depth (m),

- ``\mathcal{M}`` is the neural network,

- ``SWE`` is the Snow Water Equivalent, or the height of water if all the snow melted (m),

- ``φ`` is the relative humidity (0-1),

- ``R`` is the solar radiation (W/m²).

- ``v`` is the wind speed (W/m²).

- ``T_{air}`` is the air temperature (degrees C).

- ``P_{snow}`` is the water equivalent rate of snow precipitation (m/s).

The model is a 1D model to permit utilization over any desired grid resolution and shape.

We will use the forcings and snow depth data as a validation
of the model, so the initial conditions will be the initial value
provided in the existing data.

We begin by importing the developed code to create and run the neural network,
as well as some preliminary packages:

````julia
using DataFrames, CSV, HTTP, Dates
using ClimaLand.Snow.ModelTools
using ClimaLand.Snow.DataTools
````

and also some purpose-made functions for analyzing and displaying the output.

````julia
using ClimaLand
code_dir = joinpath(pkgdir(ClimaLand), "docs/tutorials/standalone/Snow")
include(joinpath(code_dir, "analysis_tools.jl"))
include(joinpath(code_dir, "display_tools.jl"));
````

Next, we set up values of the network hyperparameters, including the
number of epochs to train it, as well as the width parameter ``n`` as outlined
in the associated parameter, and the two loss function hyperparameters ``n_1``, ``n_2``.

````julia
n = 5
n1 = 1
n2 = 1;
````

We next outline which variables in the dataset will be used as predictors,
calling them by their column name as a `Symbol`. The number and choice of these
can be changed to reflect any dataset. Another column is specified as the target
variable, in this case, the ``\frac{dz}{dt}`` column.

````julia
pred_vars = [
    :z,
    :SWE,
    :rel_hum_avg,
    :sol_rad_avg,
    :wind_speed_avg,
    :air_temp_avg,
    :dprecipdt_snow,
]
target = :dzdt;
````

Specifying the indices of the depth and precipitation variables
(used in the constraints) and the total number of input features
will be necessary when creating the model, so we will specify them
here as well.

````julia
nfeatures = length(pred_vars)
z_idx = 1
p_idx = 7;
````

We next read in the already-cleaned dataset, though for custom datasets
there is plenty of functionality provided in [`DataTools`](@ref
ClimaLand.Snow.DataTools) to scrape SNOTEL data directly. We also set the
unit timestep seen in this data (daily, so 1 day) to be used for
setting the network's constraints as well as generating timeseries during usage.
To see the code that generated this data file, check out the [data tutorial](../data_tutorial/).

````julia
data_download_link = "https://caltech.box.com/shared/static/n59m3iqcgr60gllp65rsrd3k0mtnsfmg.csv"
data = CSV.read(HTTP.get(data_download_link).body, DataFrame)
Δt = Second(86400);
````

With this, we can begin the actual usage pipeline. First, we split the
precipitation feature into rain and snow constituents, and apply a set of physically-informed
filters before extracting the necessary features with `prep_data`:

````julia
usedata = prep_data(data);
````

After this, we determine scalings for the input and target data
that are conducive to beneficial weight updates. In this case, the
target data during training will be scaled in the -1 to 1 range, and
the neural network will scale input features according to their
standard deviations (no shifting is carried out in this case, so that
the physical meaning of "0" is preserved). This data
is then converted into matrix form for ease of its conversion into
a Flux `DataLoader` object, later, during training.

````julia
out_scale = maximum(abs.(usedata[!, target]))
in_scales = std.(eachcol(select(usedata, pred_vars)))
x_train, y_train = make_data(usedata, pred_vars, target, out_scale);
````

We then create the model itself given the hyperparameters specified
above, and indicate which features are to be used to determine the
boundary constraints on the network, and return the trainable weights
for the overall model.

````julia
model = make_model(nfeatures, n, z_idx, p_idx, in_scale = in_scales)
ps = get_model_ps(model);
````

As training updates are better with the scaled data, we have to modify
the timescale and output scaling of the model structure prior to training.
This step is undone/reset after training is over.

````julia
settimescale!(model, Dates.value(Δt) * out_scale)
setoutscale!(model, 1.0);
````

With that, training is as simple as calling the `trainmodel!` function:

````julia
print("\nTraining model!\n")
trainmodel!(model, ps, x_train, y_train, n1, n2, verbose = true);
````

````

Training model!
Epoch: 10 | training loss: 0.022190021
Epoch: 20 | training loss: 0.019470619
Epoch: 30 | training loss: 0.019247133
Epoch: 40 | training loss: 0.019255592
Epoch: 50 | training loss: 0.019042684
Epoch: 60 | training loss: 0.018979127
Epoch: 70 | training loss: 0.019085307
Epoch: 80 | training loss: 0.018973015
Epoch: 90 | training loss: 0.01938594
Epoch: 100 | training loss: 0.019005865

````

To show the model's output on some of our training data in physically meaningful
units, we first reset the timesacle and output scaling constants. From there,
all we do is pass the dataframe for a given SNOTEL site and the trained model
to the `make_timeseries` function, and we can compare the result to the actual data.

````julia
print("\nGenerating Timeseries!\n")
setoutscale!(model, out_scale)
settimescale!(model, Dates.value(Δt));
````

````

Generating Timeseries!

````

For instance, let's show the results on site 1286 (Slagamount Lakes site, Montana):

````julia
site_id = 1286
sitedata = usedata[usedata[!, :id] .== site_id, :]
true_series = sitedata[!, :z]
pred_series, _, _ = make_timeseries(model, sitedata, Δt)
display_scores(pred_series, true_series, timeseries = true)
siteplot(
    sitedata[!, :date],
    pred_series,
    true_series,
    "SITE: " * string(site_id) * ", Snow Depth (m)",
    savename = "base_tutorial_plot.png",
)
````

````
"/Users/mitraasadollahi/Projects/CLIMA_2024/ClimaLand.jl/docs/src/generated/standalone/Snow/base_tutorial_plot.png"
````

![](base_tutorial_plot.png)

Additional functionality can be explored through the optional arguments
to the developed functions, though creating timeseries for any validation
dataset can be handled with a similar call to `make_timeseries`. The timestep
`Δt` (as well as a matching call to the network with `settimescale!`) can
also be changed to different values to evaluate the network's capability on
validation data with different temporal resolutions, without the need
for retraining.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

