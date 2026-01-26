"""
    RegressionModel <: StatisticalModel

Abstract supertype for all regression models.
"""
abstract type RegressionModel <: StatisticalModel end

"""
    fitted(model::RegressionModel)

Return the fitted values of the model.
"""
function fitted end

"""
    response(model::RegressionModel)

Return the model response (a.k.a. the dependent variable).
"""
function response end

"""
    responsename(model::RegressionModel)

Return the name of the model response (a.k.a. the dependent variable).
"""
function responsename end

"""
    meanresponse(model::RegressionModel)

Return the mean of the response.
"""
function meanresponse end

"""
    modelmatrix(model::RegressionModel)

Return the model matrix (a.k.a. the design matrix).
"""
function modelmatrix end

"""
    crossmodelmatrix(model::RegressionModel)

Return `X'X` where `X` is the model matrix of `model`.
This function will return a pre-computed matrix stored in `model` if possible.
"""
crossmodelmatrix(model::RegressionModel) = (x = modelmatrix(model); Symmetric(x' * x))

"""
    leverage(model::RegressionModel)

Return the diagonal of the projection matrix of the model.
"""
function leverage end

"""
    cooksdistance(model::RegressionModel)

Compute [Cook's distance](https://en.wikipedia.org/wiki/Cook%27s_distance)
for each observation in linear model `model`, giving an estimate of the influence
of each data point.
"""
function cooksdistance end

"""
    residuals(model::RegressionModel)

Return the residuals of the model.
"""
function residuals end

"""
    predict(model::RegressionModel, [newX])

Form the predicted response of `model`. An object with new covariate values `newX` can be supplied,
which should have the same type and structure as that used to fit `model`; e.g. for a GLM
it would generally be a `DataFrame` with the same variable names as the original predictors.
"""
function predict end

"""
    predict!

In-place version of [`predict`](@ref).
"""
function predict! end

"""
    dof_residual(model::RegressionModel)

Return the residual degrees of freedom of the model.
"""
function dof_residual end

"""
    reconstruct(model::RegressionModel[, newY])

Reconstruct explanatory variables from `model`.
An object with new response values `newX` can be supplied, which should have
the same type and structure as the output of [`predict(model)`](@ref).
"""
function reconstruct end

"""
    reconstruct!

In-place version of [`reconstruct`](@ref).
"""
function reconstruct! end

"""
    offset(model::RegressionModel)

Return the offset used in the model, i.e. the term added to the linear predictor with
known coefficient 1, or `nothing` if the model was not fit with an offset.
"""
function offset end

"""
    linearpredictor(model::RegressionModel)

Return the model's linear predictor, `Xβ` where `X` is the model matrix and `β` is the
vector of coefficients, or `Xβ + offset` if the model was fit with an offset.
"""
function linearpredictor end

"""
    linearpredictor!(storage, model::RegressionModel)

In-place version of [`linearpredictor`](@ref), storing the result in `storage`.
"""
function linearpredictor! end

"""
    vif(m::RegressionModel)

Compute the variance inflation factor (VIF).

The [VIF](https://en.wikipedia.org/wiki/Variance_inflation_factor) measures
the increase in the variance of a parameter's estimate in a model with multiple parameters relative to
the variance of a parameter's estimate in a model containing only that parameter.

See also [`gvif`](@ref).

!!! warning
    This method will fail if there is (numerically) perfect multicollinearity,
    i.e. rank deficiency. In that case though, the VIF
    is not particularly informative anyway.
"""
function vif end
# This generic function is owned by StatsModels.jl, which is the sole provider
# of the default definition.

"""
    gvif(m::RegressionModel; scale=false)

Compute the generalized variance inflation factor (GVIF).

If `scale=true`, then each GVIF is scaled by the degrees of freedom
for (number of coefficients associated with) the predictor: ``GVIF^(1 / (2*df))``.

The [GVIF](https://doi.org/10.2307/2290467)
measures the increase in the variance of a (group of) parameter's estimate in a model
with multiple parameters relative to the variance of a parameter's estimate in a
model containing only that parameter. For continuous, numerical predictors, the GVIF
is the same as the VIF, but for categorical predictors, the GVIF provides a single
number for the entire group of contrast-coded coefficients associated with a categorical
predictor.

See also [`vif`](@ref).

## References

Fox, J., & Monette, G. (1992). Generalized Collinearity Diagnostics.
Journal of the American Statistical Association, 87(417), 178. doi:10.2307/2290467
"""
function gvif end
# This generic function is owned by StatsModels.jl, which is the sole provider
# of the default definition.
