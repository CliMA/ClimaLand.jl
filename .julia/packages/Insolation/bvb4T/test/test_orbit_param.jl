rtol = 1e-2
od = Insolation.OrbitalData()
@test mod(Insolation.ϖ_spline(od, 0.0), 2π) ≈
      mod(IP.lon_perihelion_epoch(param_set), 2π) rtol = rtol
@test Insolation.γ_spline(od, 0.0) ≈ IP.obliq_epoch(param_set) rtol = rtol
@test Insolation.e_spline(od, 0.0) ≈ IP.eccentricity_epoch(param_set) rtol =
    rtol

ϖ0, γ0, e0 = orbital_params(od, 0.0)
@test mod(ϖ0, 2π) ≈ mod(IP.lon_perihelion_epoch(param_set), 2π) rtol = rtol
@test γ0 ≈ IP.obliq_epoch(param_set) rtol = rtol
@test e0 ≈ IP.eccentricity_epoch(param_set) rtol = rtol
