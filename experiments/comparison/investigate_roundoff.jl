# Verify the claimed cause of the 1-ulp ρe_int difference between `Column` and
# `ColumnEnsemble` after one step: the extruded multi-column geometry carries a
# constant horizontal metric factor J_h = R²·cosd(lat)·(π/180)² in J/WJ/invJ.
# It cancels analytically inside the vertical divergence stencil
# (Ju³₊ − Ju³₋)·invJ, but each J_h·J_v product (and invJ = 1/(J_h·J_v)) rounds,
# so tendencies can land 1 ulp apart. CPU only.
#
# Run: julia --project=.buildkite experiments/comparison/investigate_roundoff.jl

# Builds soil_column / soil_ensemble, both sims, and takes ONE step of each
# (this is where the ρe_int difference is created).
include(joinpath(@__DIR__, "compare_soil.jl"))

import ClimaCore: Fields, Geometry, Operators, Spaces

# Signed difference between two Float64s measured in ulps of the larger one.
ulps(a, b) = (a - b) / eps(max(abs(a), abs(b)))

# ---------------------------------------------------------------------------
# CHECK 1: the two spaces discretize the SAME vertical mesh (z is bitwise
# identical), so any state difference must come from something else.
# ---------------------------------------------------------------------------
z_col = vec(parent(Fields.coordinate_field(soil_column.space.subsurface).z))   # 25 center z's, 1-D column
z_ens = vec(parent(Fields.coordinate_field(soil_ensemble.space.subsurface).z)) # 25 center z's, ensemble col 1
@info "CHECK 1: z bitwise identical?" all(z_col .== z_ens)                     # expect true

# ---------------------------------------------------------------------------
# CHECK 2: the Jacobians are NOT the same. The 1-D column has J = J_v (≈ Δz);
# the extruded ensemble has J = J_h · J_v with the sphere-metric factor
# J_h = R²·cosd(lat)·(π/180)² (product_geometry, globalgeometry.jl:234).
# ---------------------------------------------------------------------------
J_col = vec(parent(Fields.local_geometry_field(soil_column.space.subsurface).J))   # pure vertical J_v
J_ens =
    vec(parent(Fields.local_geometry_field(soil_ensemble.space.subsurface).J)) # J_h · J_v
J_h_expected = 6.371229e6^2 * cosd(0.0) * (π / 180)^2       # R = 6.371229e6 m (default), lat = 0 → ≈ 1.2367e10
@info "CHECK 2: J ratio (should all be ≈ J_h ≈ 1.2367e10)" extrema(
    J_ens ./ J_col,
) J_h_expected

# ---------------------------------------------------------------------------
# CHECK 3 (the mechanism, in pure Float64 arithmetic — no ClimaCore): scaling
# numerator and denominator of the divergence stencil by J_h cancels exactly
# in ℝ but not in floating point. Count how often the two evaluations differ.
# ---------------------------------------------------------------------------
J_face = vec(
    parent(Fields.local_geometry_field(soil_ensemble.space.subsurface_face).J),
) # face J's incl. J_h
n_diff = 0                                                   # how many stencil evaluations differ
max_ulp = 0.0                                                # and by how much
for i in 1:(length(J_col) - 1)                               # one divergence-like stencil per interior center
    w₊, w₋ = sin(13.0 * i), cos(17.0 * i)                    # arbitrary reproducible "flux" values
    Jv₊, Jv₋, Jvc =
        J_face[i + 1] / J_h_expected, J_face[i] / J_h_expected, J_col[i] # vertical-only Jacobians
    plain = (Jv₊ * w₊ - Jv₋ * w₋) * inv(Jvc)                 # what the 1-D Column stencil computes
    scaled =
        ((J_h_expected * Jv₊) * w₊ - (J_h_expected * Jv₋) * w₋) *
        inv(J_h_expected * Jvc)                              # what the ensemble stencil computes
    global n_diff += plain != scaled                         # bitwise different?
    global max_ulp = max(max_ulp, abs(ulps(plain, scaled)))  # size of the difference
end
@info "CHECK 3: scalar stencil with/without J_h" n_diff_of =
    (n_diff, length(J_col) - 1) max_ulp

# ---------------------------------------------------------------------------
# CHECK 4 (the mechanism, through the REAL operator): apply DivergenceF2C to
# the SAME face field on both spaces. Analytically the results are equal;
# bitwise they differ at the ~1-ulp level. THIS is where the difference
# manifests — inside vertical stencil evaluations on the ensemble space.
# ---------------------------------------------------------------------------
wf(z) = Geometry.WVector(sin(50.0 * z))                      # smooth test flux, same formula on both spaces
w_col = wf.(Fields.coordinate_field(soil_column.space.subsurface_face).z)     # column flux
w_ens = wf.(Fields.coordinate_field(soil_ensemble.space.subsurface_face).z)   # ensemble flux (bitwise-identical input, per CHECK 1)
divf2c = Operators.DivergenceF2C()                           # the operator EnergyHydrology uses for flux divergences
d_col = vec(parent(divf2c.(w_col)))                          # divergence on the 1-D column space
d_ens = vec(parent(divf2c.(w_ens)))                          # divergence on the multi-column space
@info "CHECK 4: DivergenceF2C on identical input" n_differing =
    count(d_col .!= d_ens) max_ulps = maximum(abs.(ulps.(d_col, d_ens)))

# ---------------------------------------------------------------------------
# CHECK 5 (the symptom): after the one step taken in compare_soil.jl, find
# exactly which ρe_int entries differ and confirm each is a single ulp.
# ϑ_l should have NO differing entries (same relative perturbation is < ½ ulp
# of a state of size ~0.2, so it rounds away).
# ---------------------------------------------------------------------------
ρe_col = vec(parent(sim_column._integrator.u.soil.ρe_int))   # 25 values, single column
pe = parent(sim_ensemble._integrator.u.soil.ρe_int)          # ensemble parent, columns on last dim
ρe_ens = vec(collect(selectdim(pe, ndims(pe), 1)))           # column 1 of the ensemble
idx = findall(ρe_col .!= ρe_ens)                             # which vertical levels differ
@info "CHECK 5: ρe_int entries differing after 1 step" idx
for i in idx                                                 # each should be exactly ±1 ulp
    @info @sprintf(
        "  level %2d: col = %.17e  ens = %.17e  diff = %+.1f ulp",
        i,
        ρe_col[i],
        ρe_ens[i],
        ulps(ρe_col[i], ρe_ens[i]),
    )
end
ϑ_col = vec(parent(sim_column._integrator.u.soil.ϑ_l))       # water content, single column
pϑ = parent(sim_ensemble._integrator.u.soil.ϑ_l)             # ensemble parent
ϑ_ens = vec(collect(selectdim(pϑ, ndims(pϑ), 1)))            # column 1
@info "CHECK 5b: ϑ_l entries differing (expect 0)" count(ϑ_col .!= ϑ_ens)
