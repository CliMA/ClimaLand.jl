a_root = FT(13192)
a_stem = FT(515.5605)
b_root = FT(2.1079)
b_stem = FT(0.9631)
size_reservoir_leaf_moles = FT(16766.2790)
size_reservoir_stem_moles = FT(11000.8837)
K_max_root_moles = FT(12.9216)
K_max_stem_moles = FT(3.4415)
z_leaf = FT(12) # height of leaf
z_root_depths = [FT(-1.0)] # m, rooting depth
z_bottom_stem = FT(0.0)

root_domain = RootDomain{FT}(z_root_depths, [z_bottom_stem, z_leaf])
param_set = Roots.RootsParameters{FT}(
    a_root,
    b_root,
    a_stem,
    b_stem,
    size_reservoir_stem_moles,
    size_reservoir_leaf_moles,
    K_max_root_moles,
    K_max_stem_moles,
)

function leaf_transpiration(t::ft) where {ft}
    mass_mole_water = ft(0.018)
    T = ft(0.0)
    T_0 = ft(0.01 / mass_mole_water)
    if t < ft(500)
        T = T_0
    elseif t < ft(1000)
        T = ft(10 * (T_0 / 5) * (t - 500) / 500 + T_0)
    else
        T = ft(10 * (T_0 / 5) * 500 / 500 + T_0)
    end
    return T
end

p_soil0 = [FT(-0.02)]
transpiration = PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
root_extraction = PrescribedSoilPressure{FT}((t::FT) -> p_soil0)
roots = Roots.RootsModel{FT}(;
    domain = root_domain,
    param_set = param_set,
    root_extraction = root_extraction,
    transpiration = transpiration,
)

# Set system to equilibrium state by setting LHS of both ODEs to 0


function f!(F, Y)
    T0 = 0.01 / 0.018
    flow_in_stem = sum(
        flow.(
            z_root_depths,
            z_bottom_stem,
            p_soil0,
            Y[1],
            a_root,
            b_root,
            K_max_root_moles,
        ),
    )
    flow_out_stem = flow(
        z_bottom_stem,
        z_leaf,
        Y[1],
        Y[2],
        a_stem,
        b_stem,
        K_max_stem_moles,
    )
    F[1] = flow_in_stem - T0
    F[2] = flow_out_stem - T0
end

soln = nlsolve(f!, [-1.0, -0.9])
p_stem_ini = soln.zero[1]
p_leaf_ini = soln.zero[2]

theta_stem_0 = p_to_theta(p_stem_ini)
theta_leaf_0 = p_to_theta(p_leaf_ini)
y1_0 = FT(theta_stem_0 * size_reservoir_stem_moles)
y2_0 = FT(theta_leaf_0 * size_reservoir_leaf_moles)
y0 = [y1_0, y2_0]
Y, p, coords = initialize(roots)
Y.vegetation.rwc .= y0

root_ode! = make_ode_function(roots)

t0 = FT(0);
tf = FT(60 * 60.0 * 10);
dt = FT(1);

prob = ODEProblem(root_ode!, Y, (t0, tf), p);
sol = solve(prob, Euler(), dt = dt);

dY = similar(Y)
root_ode!(dY, Y, p, 0.0)
@test sqrt(mean(dY.vegetation.rwc .^ 2.0)) < 1e-8 # starts in equilibrium


y_1 = reduce(hcat, sol.u)[1, :]
y_2 = reduce(hcat, sol.u)[2, :]
y_theta_1 = y_1 ./ size_reservoir_stem_moles
y_theta_2 = y_2 ./ size_reservoir_leaf_moles
p_stem = theta_to_p.(y_theta_1)
p_leaf = theta_to_p.(y_theta_2)

function f2!(F, Y)
    p_soilf = p_soil0
    Tf = 0.01 / 0.018 .* 3.0
    flow_in_stem = sum(
        flow.(
            z_root_depths,
            z_bottom_stem,
            p_soilf,
            Y[1],
            a_root,
            b_root,
            K_max_root_moles,
        ),
    )
    flow_out_stem = flow(
        z_bottom_stem,
        z_leaf,
        Y[1],
        Y[2],
        a_stem,
        b_stem,
        K_max_stem_moles,
    )
    F[1] = flow_in_stem - Tf
    F[2] = flow_out_stem - Tf
end


# Check that the final state is in the new equilibrium
soln = nlsolve(f2!, [-1.0, -0.9]; ftol = 1e-10)
p_stem_f = soln.zero[1]
p_leaf_f = soln.zero[2]
@test abs(p_stem_f - p_stem[end]) < 1e-10
@test abs(p_leaf_f - p_leaf[end]) < 1e-10


#= Plots for comparison to Anna's script. Can remove when ready to merge.
using Plots
times = sol.t .<=60*60.0*2
plot(sol.t[times], y_1[times], label="stem", xaxis="t [s]", yaxis="water content [mol]", dpi=500)
plot!(sol.t[times], y_2[times], label="leaf")
plot!(sol.t[times], y1_0.*ones(sum(times),1),label="W_stem_0",dpi=500)
plot!(sol.t[times], y2_0.*ones(sum(times),1),label="W_leaf_0",dpi=500)
savefig("absolute_water_content.png") 

plot(sol.t[times], y_theta_1[times], label="stem", xaxis="t [s]", yaxis="relative water content [m3/m3]",dpi=500)
plot!(sol.t[times], y_theta_2[times], label="leaf", dpi=500)
savefig("relative_water_content.png") 

plot(sol.t[times],p_stem[times],linewidth=2,xaxis="time [s]",yaxis="pressure [MPa]",label="stem",dpi=500)
plot!(sol.t[times],p_leaf[times],linewidth=2,label="leaf",dpi=500)
plot!(sol.t[times],p_stem_ini.*ones(sum(times),1),label="p_stem_0",dpi=500)
plot!(sol.t[times],p_leaf_ini.*ones(sum(times),1),label="p_leaf_0",dpi=500)
savefig("pressure.png") 

=#
