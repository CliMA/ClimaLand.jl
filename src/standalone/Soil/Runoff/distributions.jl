using Statistics
using StatsBase
using SpecialFunctions
using NLsolve

struct FrechetDistribution end
function constraint!(F, α, x)
    F[1] = 1/α[1]+ mean(x.^(-α[1]) .*log.(x))/mean(x.^(-α[1])) - mean(log.(x))
end


function fit(dist::FrechetDistribution, x)
    wrapper!(F, α) = constraint!(F, α, x)
    α = nlsolve(wrapper!, [1.0])
    if α.f_converged && α.zero[1] >0
        σ = mean(x.^(-α.zero[1]))^(-1/α.zero[1])
        return Statistics.mean(x), Statistics.var(x), [α.zero[1], σ]
    else
        return Statistics.mean(x), Statistics.var(x), [NaN, NaN]
    end
    
end

function pdf(dist::FrechetDistribution, x, params)
    (α, σ) = params
    y = x/σ
    return α/σ*exp(-y^(-α))*y^(-α-1)
end


n_params(dist::FrechetDistribution) = 2

struct LogNormalDistribution end

function fit(dist::LogNormalDistribution, y)
    μ = mean(log.(y))
    σ2 = mean((log.(y) .- μ).^2)
    return Statistics.mean(y), Statistics.var(y), [μ, σ2]
end

function pdf(dist::LogNormalDistribution, x, params)
    (μ, σ2) = params
    return 1/(x*sqrt(σ2*2π))*exp(-(log(x)-μ)^2/2/σ2)
end


n_params(dist::LogNormalDistribution) = 2


struct PercentileDistribution end

function fit(dist::PercentileDistribution, x; q = [0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1.0])
    x̄ = Statistics.mean(x)
    var = Statistics.var(x)
    params = StatsBase.quantile(x, q)
    return x̄, var, params
end

n_params(dist::PercentileDistribution) = 11
function sortedIndex(array, value)
    low = 0
    high = length(array)
    while low < high
        idx = low + high
        if array[idx] < value
            low = idx + 1
        else
            high -= 1
        end
    end
    return low
end

function cdf(dist::PercentileDistribution, x, p; y = Array(0.0:0.1:1.0))
    id_u = sortedIndex(p, x)
    if id_u <= 1
        return  0
    elseif id_u == length(y)+1
        return 1
    else
        return y[id_u-1] + 0.1/(p[id_u]-p[id_u-1])*(x - p[id_u-1])
    end
end

function pdf(dist::PercentileDistribution, x, p; y = Array(0.0:0.1:1.0))
    id_u = sortedIndex(p, x)
    if id_u <= 1
        return  eps(Float64)
    elseif id_u == length(y)+1
        return eps(Float64)
    else
        return 0.1/(p[id_u]-p[id_u-1])
    end
end




struct GammaDistribution end

function update_k(k, s)
    digamma_deriv = (digamma(k + 0.001 * k) - digamma(k)) / (0.001 * k)
    kguess = k - (log(k) - digamma(k) - s) / (1 / k - digamma_deriv)
    dk = kguess - k
    return kguess, dk
end

function fit(dist::GammaDistribution, x)
    x̄ = Statistics.mean(x)
    var = Statistics.var(x)
    s = log(x̄) - Statistics.mean(log.(x))
    k = (3 - s + sqrt((s - 3)^2 + 24 * s)) / (12 * s)
    dk = Inf
    while abs(dk / k) > 1e-3
        k, dk = update_k(k, s)
    end
    k += dk
    θ = x̄ / k
    β = 1 / θ
    α = k
    return x̄, var, [α, β]
end

n_params(dist::GammaDistribution) = 2
function pdf(dist::GammaDistribution, x, params)
    (α, β) = params
    return max(eps(typeof(x)), β^α / gamma(α) * (x)^(α - 1) * exp(-β * x))
end



struct InvGammaDistribution end

function fit(dist::InvGammaDistribution, y)
    # if y ~ InvGamma(α, β), 1/y ~ Gamma(α,β)
    x = 1 ./ y
    x̄ = Statistics.mean(x)
    s = log(x̄) - Statistics.mean(log.(x))
    k = (3 - s + sqrt((s - 3)^2 + 24 * s)) / (12 * s)
    dk = Inf
    while abs(dk / k) > 1e-3
        k, dk = update_k(k, s)
    end
    k += dk
    θ = x̄ / k
    β = 1 / θ
    α = k
    return Statistics.mean(y), Statistics.var(y), [α, β]
end

function pdf(dist::InvGammaDistribution, x, params)
    (α, β) = params
    return max(eps(typeof(x)), β^α / gamma(α) * (1 / x)^(α + 1) * exp(-β / x))

end

n_params(dist::InvGammaDistribution) = 2


function log_likelihood(dist, data::Vector, params)
    return sum(log.(pdf.(Ref(dist), data, Ref(params))))
end
#=
struct GeneralizedGammaDistribution end
n_params(dist::GeneralizedGammaDistribution) = 3
function beta(α, δ, x)
    return (Statistics.mean(x.^δ)/α)^(1/δ)
end

function f_of_alpha(β, δ, x)
    z = @. log(x/β)^δ
    f = δ*(Statistics.mean(log.(x)) - Statistics.mean(z .^δ .* log.(x))/Statistics.mean(x .^δ))
    return 1/f - α
end

h_of_delta(α, δ,x) = -digamma(α) + δ*Statistics.mean(log.(x)) - log(Statistics.mean(x .^δ)) + log(α)
#dhda(α, δ,x) = -(digamma(α+0.001*α)-digamma(α))/(0.001*α)+ 1/α
#dhdd(α, δ,x) = Statistics.mean(log.(x))  - 1/Statistics.mean(x .^δ)*Statistics.mean(x .^(δ-1))*δ

function fit(dist::GeneralizedGammaDistribution, x)
    x̄ = Statistics.mean(x)
    var = Statistics.var(x)
    function f!(F, Y; x = x)
        α, δ = Y
        @assert α > 0
        @assert δ > 0 
        β = beta(α, δ, x)
        @assert β > 0 
        F[1] = f_of_alpha(β, δ, x)
        F[2] = h_of_delta(α, δ, x)
    end
    Y0 = [2.0, 1.0]
    nlsolve(f!, Y0)
    return x̄, var, [α, β, δ]
end

function log_likelihood_ggd(data::Vector, parameters::Vector{FT}) where {FT}
    α, δ = parameters
    c1 = α > 0
    c2 = δ > 0
    β = beta(α, δ, data)
    if !c1 || !c2
        return -FT(1.0)/eps(FT)
    else
        return (δ*α-1)*Statistics.mean(log.(x)) - Statistics.mean((x./β).^δ)
    end
    nothing
end
function fit_ggd(data::Vector, initial_guess::Vector{FT}) where{FT}
    function wrapper(guess::Vector{FT}) where{FT}
        return -(log_likelihood_ggd(data,guess))
    end
    results = optimize(wrapper, initial_guess, NelderMead())
    return Optim.minimizer(results), Optim.minimum(results)
end

=#
