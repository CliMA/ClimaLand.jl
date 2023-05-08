struct ACM{FT <: AbstractFloat}
    d1::FT
    θ::FT
    k::FT 
    d2::FT
    b2::FT
    c1::FT
    a2::FT
    c2::FT
    eb1T::FT
    ψ_d::FT
    H::FT
end

function ACM(::Type{FT}) where {FT <: AbstractFloat}
    d1 = FT(0.0156935) 
    θ = FT(4.22273)
    k = FT(208.868)
    d2 = FT(0.0453194)
    b2 = FT(0.37836)
    c1 = FT(7.19298)
    a2 = FT(0.011136)
    c2 = FT(2.1001)
    eb1T = FT(0.789798)
    ψ_d = FT(-2)
    H = FT(1)
    return ACM(d1, θ, k, d2, b2, c1, a2, c2, eb1T, ψ_d, H)
end

function (m::ACM{FT})(; lat, doy, t_max, t_min, lai, rad, ca, ce) where {FT}

    # these parameters are directly from the ACM C files,
    # which differ from the original
    # implementation in Williams et al. (1997)
    
    # compute daily canopy conductance, gc
    gc = (abs(m.ψ_d)^m.eb1T) / (m.b2 * m.H + FT(0.5) * (t_max - t_min))
    
    # compute p parameter needed for ci
    p = lai .* FT(1) .* ce .* exp.(m.a2 .* t_max) ./ gc
 
    
    # compute the q parametefr needed for ci
    q = m.θ .- m.k
    
    # compute the internal CO2 concentration, ci
    ci = FT(0.5) .* (ca .+ q .- p .+ sqrt.((ca .+ q .- p).^FT(2) .- FT(4) .* (ca .* q .- p .* m.θ)))
    #println("ACM ci: ", ci)
    
    # compute canopy-level quantum yield, e0
    e0 = m.c1 .* (lai .^ FT(2)) ./ (m.c2 .+ lai .^ FT(2))
    
    # compute the day length dayl
    dec = FT(-23.4) * cos((FT(360) * (doy + FT(10)) / FT(365)) * FT(π) / FT(180)) * FT(π) / FT(180)
    mult = tan(lat * FT(π) / FT(180)) * tan(dec)
    if mult >= FT(1)
        dayl = FT(24)  
    elseif mult <= -1 
        dayl = FT(0)
    else
      dayl = FT(24) * FT(acos(-mult)) / FT(π)
    end
    
    # compute co2 rate of diffusion to the site of fixation, pd
    pd = gc .* (ca .- ci)
    
    # compute light limitation pi
    pi = e0 .* rad .* pd ./ (e0 .* rad .+ pd)
    
    # compute gpp
    gpp = pi .* (m.d1 .* dayl .+ m.d2)    
    return gpp
end