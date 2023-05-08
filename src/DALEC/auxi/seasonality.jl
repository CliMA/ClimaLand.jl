"""

lab_release_factor

    Compute the labile release factor.

## Arguments:
- `t::Real`: time (in days)
- `lab_lifespan`: labile lifespan
- `clab_release_period`: Clab Release period
- `Bday`: Bday
- `FT::Type`: Float type, either Float32 or Float64
"""
function lab_release_factor(t, lab_lifespan, clab_release_period, Bday, FT)
    fl = (log(lab_lifespan) - log(lab_lifespan - 1)) * FT(0.5)
    wl = clab_release_period * sqrt(FT(2)) / FT(2)
    osl = offset(lab_lifespan, wl, FT)
    sf = FT(365.25) / FT(π)
    return (FT(2) / sqrt(FT(π))) * (fl / wl) * exp( - (sin((t - Bday + osl) / sf) * sf / wl)^2)
end;


"""

leaf_fall_factor

    Compute the leaf fall factor.

## Arguments:
- `t::Real`: time (in days)
- `leaf_lifespan`: leaf lifespan
- `leaf_fall_period`: leaf fall period
- `Fday`: Fday
- `FT::Type`: Float type, either Float32 or Float64
"""
function leaf_fall_factor(t, leaf_lifespan, leaf_fall_period, Fday, FT)
    ff = (log(leaf_lifespan) - log(leaf_lifespan - 1)) * FT(0.5)
    wf = leaf_fall_period * sqrt(FT(2)) / 2

    osf = offset(leaf_lifespan, wf, FT)
    sf = FT(365.25)/FT(π)
    return (2/sqrt(FT(π)))*(ff/wf)*exp(-(sin((t-Fday+osf)/sf)*sf/wf)^2)
end

function offset(L, w, FT)
    mxc = [FT(0.000023599784710),FT(0.000332730053021),FT(0.000901865258885),FT(-0.005437736864888),FT(-0.020836027517787), FT(0.126972018064287),FT(-0.188459767342504)]
    lf = log(L-1)
    os = mxc[1]*lf^6+mxc[2]*lf^5+mxc[3]*lf^4+mxc[4]*lf^3+mxc[5]*lf^2+mxc[6]*lf+mxc[7]
    os = os * w
    return os
end;