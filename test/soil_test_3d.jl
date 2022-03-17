ν = FT(0.495);
Ksat = FT(1)#0.0443 / 3600 / 100); # m/s
S_s = FT(1)#e-3); #inverse meters
vg_n = FT(2.0);
vg_α = FT(2.6); # inverse meters
vg_m = FT(1) - FT(1) / vg_n;
θ_r = FT(0);
zmax = FT(0);
zmin = FT(-1);
xmax = FT(1.0)
soil_domain = HybridBox(
    Float64;
    xlim = (0.0, xmax),
    ylim = (0.0, 1.0),
    zlim = (zmin, zmax),
    nelements = (500, 2, 500),
    npolynomial = 3,
)
top_flux_bc = FT(0.0);
bot_flux_bc = FT(0.0);
sources = ()
boundary_fluxes = FluxBC{FT}(top_flux_bc, bot_flux_bc)
params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r);

soil = Soil.RichardsModel{FT}(;
    param_set = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = sources,
)

# sinusoidal water table
function init_soil!(Ysoil, x, z, params)
    function hydrostatic_profile(
        x::FT,
        z::FT,
        params::RichardsParameters{FT},
    ) where {FT}
        @unpack ν, vg_α, vg_n, vg_m, θ_r, S_s = params
        z_∇ = FT(zmin / 2.0 + (zmax - zmin) / 10.0 * sin(π * 2 * x / xmax))
        if z > z_∇
            S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
            ϑ_l = S * (ν - θ_r) + θ_r
        else
            ϑ_l = -S_s * (z - z_∇) + (ν - θ_r)
        end
        return FT(ϑ_l)
    end
    Ysoil.soil.ϑ_l .= hydrostatic_profile.(x, z, Ref(params))
end
soil_ode! = make_ode_function(soil)
Y, p, coords = initialize(soil)
init_soil!(Y, coords.x, coords.z, soil.param_set)
dY = similar(Y)
soil_ode!(dY, Y, p, 0.0)

function dθdx(x, z)
    z_∇ = zmin / 2.0 + (zmax - zmin) / 10.0 * sin(π * 2 * x / xmax)
    dz∇dx = (zmax - zmin) / 10.0 * cos(π * 2 * x / xmax) * (2 * π / xmax)
    if z > z_∇
        dSdz_∇ =
            -(vg_m) *
            (1 + (vg_α * (z - z_∇))^vg_n)^(-vg_m - 1) *
            vg_n *
            (vg_α * (z - z_∇))^(vg_n - 1) *
            (-vg_α)
        dSdx = dSdz_∇ * dz∇dx
        dθdx = (ν - θ_r) * dSdx
    else
        dθdx = S_s * dz∇dx
    end
    return dθdx
end

function d2θdx2(x, z)

    z_∇ = zmin / 2.0 + (zmax - zmin) / 10.0 * sin(π * 2 * x / xmax)
    dz∇dx = (zmax - zmin) / 10.0 * cos(π * 2 * x / xmax) * (2π / xmax)
    d2z∇dx2 = -(zmax - zmin) / 10.0 * sin(π * 2 * x / xmax) * (2π / xmax)^2.0
    if z > z_∇
        dSdz_∇ =
            -(vg_m) *
            (1 + (vg_α * (z - z_∇))^vg_n)^(-vg_m - 1) *
            vg_n *
            (vg_α * (z - z_∇))^(vg_n - 1) *
            (-vg_α)

        d2Sdz∇2 =
            -vg_m *
            vg_n *
            vg_α^2.0 *
            (
                (vg_n - 1) *
                (vg_α * (z - z_∇))^(vg_n - 2) *
                (1 + (vg_α * (z - z_∇))^vg_n)^(-vg_m - 1) -
                (vg_m + 1) *
                vg_n *
                (vg_α * (z - z_∇))^(2.0 * (vg_n - 1)) *
                (1 + (vg_α * (z - z_∇))^vg_n)^(-vg_m - 2)
            )
        d2Sdx2 = dz∇dx^2.0 * d2Sdz∇2 + dSdz_∇ * d2z∇dx2
        d2θdx2 = (ν - θ_r) * d2Sdx2
    else
        d2θdx2 = S_s * d2z∇dx2
    end
    return d2θdx2
end

function dKdθ(θ)
    S = (θ - θ_r) / (ν - θ_r)
    if S < 1
        f = 1.0 - (1.0 - S^(1.0 / vg_m))^vg_m
        f1 = f^2.0 / 2.0 / S^0.5
        f2 = 2 * S^(1 / vg_m - 1 / 2) * f / (1 - S^(1 / vg_m))^(1.0 - vg_m)
        return (f1 + f2) * Ksat / (ν - θ_r)
    else
        return 0.0
    end

end

function dψdθ(θ)
    S = (θ - θ_r) / (ν - θ_r)
    if S < 1.0
        return 1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r) *
               (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
               S^(-1 / vg_m - 1)
    else
        return 1.0 / S_s
    end
end


function d2ψdθ2(θ)
    S = (θ - θ_r) / (ν - θ_r)
    if S < 1.0
        return 1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r)^2.0 * (
            S^(-2.0 / vg_m - 2.0) *
            (-1 / vg_m) *
            (1 / vg_n - 1) *
            (S^(-1 / vg_m) - 1)^(1 / vg_n - 2) +
            (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
            (-1 / vg_m - 1) *
            S^(-1 / vg_m - 2)
        )
    else
        return 0.0
    end
end

function K(θ)
    S = (θ - θ_r) / (ν - θ_r)
    if S < 1
        f = 1.0 - (1.0 - S^(1.0 / vg_m))^vg_m
        return Ksat * f^2.0 * sqrt(S)
    else
        return Ksat
    end

end


X = coords.x
Z = coords.z
θ = Y.soil.ϑ_l

#5.517176201418359e-9 max error with horizontal terms off

# Test in unsaturated zone
myXslice = parent(X)[350, 1, 1, 1, :]
myθslice = parent(Y.soil.ϑ_l)[350, 1, 1, 1, :]
mydYslice = parent(dY.soil.ϑ_l)[350, 1, 1, 1, :]
p = sortperm(myXslice)

Xsort = myXslice[p]
θsort = myθslice[p]
dYsort = mydYslice[p]
unique_indices = unique(i -> Xsort[i], 1:length(Xsort))
XX = Xsort[unique_indices]
θX = θsort[unique_indices]
myZ = unique(parent(Z)[350, 1, 1, 1, :])[1]
dYX = dYsort[unique_indices]
expected = (@. dKdθ(θX) * dψdθ(θX) * dθdx(XX, myZ)^2.0 +
    K(θX) * dψdθ(θX) * d2θdx2(XX, myZ) +
    K(θX) * dθdx(XX, myZ)^2.0 * d2ψdθ2(θX))
@test maximum(abs.(dYX .- expected)) < 1e-8

# test in saturated zone
myXslice = parent(X)[150, 1, 1, 1, :]
myθslice = parent(Y.soil.ϑ_l)[150, 1, 1, 1, :]
mydYslice = parent(dY.soil.ϑ_l)[150, 1, 1, 1, :]
p = sortperm(myXslice)

Xsort = myXslice[p]
θsort = myθslice[p]
dYsort = mydYslice[p]
unique_indices = unique(i -> Xsort[i], 1:length(Xsort))
XX = Xsort[unique_indices]
θX = θsort[unique_indices]
myZ = unique(parent(Z)[150, 1, 1, 1, :])[1]
dYX = dYsort[unique_indices]
expected = (@. dKdθ(θX) * dψdθ(θX) * dθdx(XX, myZ)^2.0 +
    K(θX) * dψdθ(θX) * d2θdx2(XX, myZ) +
    K(θX) * dθdx(XX, myZ)^2.0 * d2ψdθ2(θX))
@test maximum(abs.(dYX .- expected)) < 1e-8
