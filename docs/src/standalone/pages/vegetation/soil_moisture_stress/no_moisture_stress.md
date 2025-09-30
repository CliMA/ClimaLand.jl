# no moisture stress

While the [Piecewise  model](@ref "Piecewise Soil Moisture Stress") and
[Tuzet model](@ref "Tuzet Soil Moisture (Water Potential) Stress") have a scaler, $\beta$,
that multiply leaf net assimilation and stomatal conductance by 0 to 1 depending on moisture stress,
we also allow the option to use no sensitivity to moisture stress ($\beta$ is always equal to 1).
