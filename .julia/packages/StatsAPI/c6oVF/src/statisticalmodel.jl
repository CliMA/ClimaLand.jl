"""
    StatisticalModel

Abstract supertype for all statistical models.
"""
abstract type StatisticalModel end

"""
    coef(model::StatisticalModel)

Return the coefficients of the model.
"""
function coef end

"""
    coefnames(model::StatisticalModel)

Return the names of the coefficients.
"""
function coefnames end

"""
    coeftable(model::StatisticalModel; level::Real=0.95)

Return a table with coefficients and related statistics of the model.
`level` determines the level for confidence intervals (by default, 95%).

The returned `CoefTable` object implements the
[Tables.jl](https://github.com/JuliaData/Tables.jl/) interface, and can be
converted e.g. to a `DataFrame` via `using DataFrames; DataFrame(coeftable(model))`.
"""
function coeftable end

"""
    confint(model::StatisticalModel; level::Real=0.95)

Compute confidence intervals for coefficients, with confidence level `level` (by default 95%).
"""
function confint end

"""
    deviance(model::StatisticalModel)

Return the deviance of the model relative to a reference, which is usually when applicable
the saturated model. It is equal, *up to a constant*, to ``-2 \\log L``, with ``L``
the likelihood of the model.
"""
function deviance end

"""
    islinear(model::StatisticalModel)

Indicate whether the model is linear.
"""
function islinear end

"""
    nulldeviance(model::StatisticalModel)

Return the deviance of the null model, obtained by dropping all
independent variables present in `model`.

If `model` includes an intercept, the null model is the one with only the intercept;
otherwise, it is the one without any predictor (not even the intercept).
"""
function nulldeviance end

"""
    loglikelihood(model::StatisticalModel)
    loglikelihood(model::StatisticalModel, observation)

Return the log-likelihood of the model.

With an `observation` argument, return the contribution of `observation` to the
log-likelihood of `model`.

If `observation` is a `Colon`, return a vector of each observation's contribution
to the log-likelihood of the model. In other words, this is the vector of the
pointwise log-likelihood contributions.

In general, `sum(loglikehood(model, :)) == loglikelihood(model)`.
"""
function loglikelihood end

"""
    nullloglikelihood(model::StatisticalModel)

Return the log-likelihood of the null model, obtained by dropping all
independent variables present in `model`.

If `model` includes an intercept, the null model is the one with only the intercept;
otherwise, it is the one without any predictor (not even the intercept).
"""
function nullloglikelihood end

"""
    score(model::StatisticalModel)

Return the score of the model, that is the gradient of the
log-likelihood with respect to the coefficients.
"""
function score end

"""
    nobs(model::StatisticalModel)

Return the number of independent observations on which the model was fitted. Be careful
when using this information, as the definition of an independent observation may vary
depending on the model, on the format used to pass the data, on the sampling plan
(if specified), etc.
"""
function nobs end

"""
    dof(model::StatisticalModel)

Return the number of degrees of freedom consumed in the model, including
when applicable the intercept and the distribution's dispersion parameter.
"""
function dof end

"""
    mss(model::StatisticalModel)

Return the model sum of squares.
"""
function mss end

"""
    rss(model::StatisticalModel)

Return the residual sum of squares of the model.
"""
function rss end

"""
    informationmatrix(model::StatisticalModel; expected::Bool = true)

Return the information matrix of the model. By default the Fisher information matrix
is returned, while the observed information matrix can be requested with `expected = false`.
"""
function informationmatrix end

"""
    stderror(model::StatisticalModel)

Return the standard errors for the coefficients of the model.
"""
stderror(model::StatisticalModel) = sqrt.(diag(vcov(model)))

"""
    vcov(model::StatisticalModel)

Return the variance-covariance matrix for the coefficients of the model.
"""
function vcov end

"""
    weights(model::StatisticalModel)

Return the weights used in the model.
"""
function weights end

"""
    isfitted(model::StatisticalModel)

Indicate whether the model has been fitted.
"""
function isfitted end

"""
Fit a statistical model.
"""
function fit end

"""
Fit a statistical model in-place.
"""
function fit! end

"""
    aic(model::StatisticalModel)

Akaike's Information Criterion, defined as ``-2 \\log L + 2k``, with ``L`` the likelihood
of the model, and `k` its number of consumed degrees of freedom
(as returned by [`dof`](@ref)).
"""
aic(model::StatisticalModel) = -2loglikelihood(model) + 2dof(model)

"""
    aicc(model::StatisticalModel)

Corrected Akaike's Information Criterion for small sample sizes (Hurvich and Tsai 1989),
defined as ``-2 \\log L + 2k + 2k(k-1)/(n-k-1)``, with ``L`` the likelihood of the model,
``k`` its number of consumed degrees of freedom (as returned by [`dof`](@ref)),
and ``n`` the number of observations (as returned by [`nobs`](@ref)).
"""
function aicc(model::StatisticalModel)
    k = dof(model)
    n = nobs(model)
    -2loglikelihood(model) + 2k + 2k*(k+1)/(n-k-1)
end

"""
    bic(model::StatisticalModel)

Bayesian Information Criterion, defined as ``-2 \\log L + k \\log n``, with ``L``
the likelihood of the model,  ``k`` its number of consumed degrees of freedom
(as returned by [`dof`](@ref)), and ``n`` the number of observations
(as returned by [`nobs`](@ref)).
"""
bic(model::StatisticalModel) = -2loglikelihood(model) + dof(model)*log(nobs(model))

function r2 end

"""
    r2(model::StatisticalModel)
    r²(model::StatisticalModel)

Coefficient of determination (R-squared).

For a linear model, the R² is defined as ``ESS/TSS``, with ``ESS`` the explained sum of squares
and ``TSS`` the total sum of squares.
"""
r2(model::StatisticalModel)

"""
    r2(model::StatisticalModel, variant::Symbol)
    r²(model::StatisticalModel, variant::Symbol)

Pseudo-coefficient of determination (pseudo R-squared).

For nonlinear models, one of several pseudo R² definitions must be chosen via `variant`.
Supported variants are:
- `:McFadden` (a.k.a. likelihood ratio index), defined as ``1 - \\log (L)/\\log (L_0)``;
- `:CoxSnell`, defined as ``1 - (L_0/L)^{2/n}``;
- `:Nagelkerke`, defined as ``(1 - (L_0/L)^{2/n})/(1 - L_0^{2/n})``.
- `:devianceratio`, defined as ``1 - D/D_0``.

In the above formulas, ``L`` is the likelihood of the model,
``L_0`` is the likelihood of the null model (the model with only an intercept),
``D`` is the deviance of the model (from the saturated model),
``D_0`` is the deviance of the null model,
``n`` is the number of observations (given by [`nobs`](@ref)).

The Cox-Snell and the deviance ratio variants both match the classical definition of R²
for linear models.
"""
function r2(model::StatisticalModel, variant::Symbol)
    loglikbased = (:McFadden, :CoxSnell, :Nagelkerke)
    if variant in loglikbased
        ll = loglikelihood(model)
        ll0 = nullloglikelihood(model)
        if variant == :McFadden
            1 - ll/ll0
        elseif variant == :CoxSnell
            1 - exp(2 * (ll0 - ll) / nobs(model))
        elseif variant == :Nagelkerke
            (1 - exp(2 * (ll0 - ll) / nobs(model))) / (1 - exp(2 * ll0 / nobs(model)))
        end
    elseif variant == :devianceratio
        dev  = deviance(model)
        dev0 = nulldeviance(model)
        1 - dev/dev0
    else
        throw(ArgumentError("variant must be one of $(join(loglikbased, ", ")) or :devianceratio"))
    end
end

const r² = r2

function adjr2 end

"""
    adjr2(model::StatisticalModel)
    adjr²(model::StatisticalModel)

Adjusted coefficient of determination (adjusted R-squared).

For linear models, the adjusted R² is defined as ``1 - (1 - (1-R^2)(n-1)/(n-p))``, with ``R^2``
the coefficient of determination, ``n`` the number of observations, and ``p`` the number of
coefficients (including the intercept). This definition is generally known as the Wherry Formula I.
"""
adjr2(model::StatisticalModel)

"""
    adjr2(model::StatisticalModel, variant::Symbol)
    adjr²(model::StatisticalModel, variant::Symbol)

Adjusted pseudo-coefficient of determination (adjusted pseudo R-squared).
For nonlinear models, one of the several pseudo R² definitions must be chosen via `variant`.
The only currently supported variants are `:McFadden`, defined as ``1 - (\\log (L) - k)/\\log (L0)`` and
`:devianceratio`, defined as ``1 - (D/(n-k))/(D_0/(n-1))``.
In these formulas, ``L`` is the likelihood of the model, ``L0`` that of the null model
(the model including only the intercept), ``D`` is the deviance of the model,
``D_0`` is the deviance of the null model, ``n`` is the number of observations (given by [`nobs`](@ref)) and
``k`` is the number of consumed degrees of freedom of the model (as returned by [`dof`](@ref)).
"""
function adjr2(model::StatisticalModel, variant::Symbol)
    k = dof(model)
    if variant == :McFadden
        ll = loglikelihood(model)
        ll0 = nullloglikelihood(model)
        1 - (ll - k)/ll0
    elseif variant == :devianceratio
        n = nobs(model)
        dev  = deviance(model)
        dev0 = nulldeviance(model)
        1 - (dev*(n-1))/(dev0*(n-k))
    else
        throw(ArgumentError("variant must be one of :McFadden or :devianceratio"))
    end
end

const adjr² = adjr2
