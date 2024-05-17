# Canopy Geometry

"""
    clumping_factor!(can::Canopy4RT{FT}, angles::SolarAngles{FT}) where {FT<:AbstractFloat}

Calculate the clumping factor, given
- `can` [`Canopy4RT`](@ref) type struct
- `angles` [`SolarAngles`](@ref) type struct
"""
function clumping_factor!(can::Canopy4RT{FT}, angles::SolarAngles{FT}) where {FT<:AbstractFloat}
    (; clump_a, clump_b) = can;
    (; sza) = angles;

    if clump_b > 0
        can.Ω = clump_a + clump_b * (1 - cosd(sza));
    end

    return nothing
end

"""
    canopy_geometry!(can::Canopy4RT{FT}, angles::SolarAngles{FT}, can_opt::CanopyOpticals{FT}, rt_con::RTCache{FT}) where {FT<:AbstractFloat}

Computes canopy optical properties (extinction coefficients for direct and
    diffuse light) based on the SAIL model. Most important input parameters are
    leaf inclination and azimuth distribution functions and sun-sensor
    geometry. Canopy clumping Ω is implemented as in Pinty et al (2015), given
- `can` [`Canopy4RT`](@ref) type struct
- `angles` [`SolarAngles`](@ref) type struct
- `can_opt` [`CanopyOpticals`](@ref) type struct
- `rt_con` [`RTCache`](@ref) type cache
"""
function canopy_geometry!(can::Canopy4RT{FT}, angles::SolarAngles{FT}, can_opt::CanopyOpticals{FT}, rt_con::RTCache{FT}) where {FT<:AbstractFloat}
    # 1. update clumping factor from zenith angle
    clumping_factor!(can, angles);

    # 2. update solor angle dependent variables
    (; sza, vza, raa) = angles;
    cos_vza = cosd(vza);
    tan_vza = tand(vza);
    cos_raa = cosd(raa);
    sin_vza = sind(vza);
    cos_sza = cosd(sza);
    sin_sza = sind(sza);
    tan_sza = tand(sza);
    cts_cto = cos_sza * cos_vza;
    dso     = sqrt( tan_sza^2 + tan_vza^2 - 2*tan_sza*tan_vza*cos_raa );
    psi_vol = abs( raa - 360*round(raa/360) );

    # 3. unpack canopy parameters
    (; dx, hot, LAI, lazitab, lidf, litab, xl, Ω) = can;

    # 4. update the RTCache
    can.cos_philo .= cosd.(lazitab .- raa);
    (; cos_philo, cos_ttli, cos_ttlo, sin_ttli) = can;

    # 5. calculate geometric factors associated with extinction and scattering
    can_opt.ks  = 0;
    can_opt.ko  = 0;
    can_opt.bf  = 0;
    can_opt.sob = 0;
    can_opt.sof = 0;
    @inbounds for i in eachindex(litab)
        _lit = litab[i];
        _lid = lidf[i];
        _ctl = cos_ttli[i];

        # interception parameters and ref/trans multipliers
        volscatt!(can.vol_scatt, sza, vza, psi_vol, _lit);
        chi_s, chi_o, frho, ftau = can.vol_scatt;

        # Extinction coefficients
        ksli = abs(chi_s / cos_sza);
        koli = abs(chi_o / cos_vza);

        # Area scattering coefficient fractions
        if ftau < 0
            sobli = abs(ftau) * FT(pi) / cts_cto;
            sofli = abs(frho) * FT(pi) / cts_cto;
        else
            sobli = frho * FT(pi) / cts_cto;
            sofli = ftau * FT(pi) / cts_cto;
        end

        # add up the values in each layer
        can_opt.ks  += ksli   * _lid;
        can_opt.ko  += koli   * _lid;
        can_opt.bf  += _ctl^2 * _lid;
        can_opt.sob += sobli  * _lid;
        can_opt.sof += sofli  * _lid;
    end

    # 6. geometric factors to be used later with rho and tau
    (; bf, ko, ks) = can_opt;
    can_opt.sdb = (ks + bf) / 2;
    can_opt.sdf = (ks - bf) / 2;
    can_opt.dob = (ko + bf) / 2;
    can_opt.dof = (ko - bf) / 2;
    can_opt.ddb = (1  + bf) / 2;
    can_opt.ddf = (1  - bf) / 2;

    # 7. eq 19 in vdT 2009 page 305 modified by Joris
    cg_con = rt_con.cg_con;
    cg_con._Cs .= cos_ttli .* cos_sza; # [nli]
    cg_con._Ss .= sin_ttli .* sin_sza; # [nli]
    cg_con._Co .= cos_ttli .* cos_vza; # [nli]
    cg_con._So .= sin_ttli .* sin_vza; # [nli]
    (; _Co, _Cs, _So, _Ss, _1s) = cg_con;
    # cg_con._cds .= _Cs * _1s .+ _Ss * cos_ttlo' ; # [nli, nlazi]
    # cg_con._cdo .= _Co * _1s .+ _So * cos_philo'; # [nli, nlazi]
    mul!(cg_con._cds, _Cs, _1s       );
    mul!(cg_con._2d , _Ss, cos_ttlo' );
    cg_con._cds .+= cg_con._2d;
    mul!(cg_con._cdo, _Co, _1s       );
    mul!(cg_con._2d , _So, cos_philo');
    cg_con._cdo .+= cg_con._2d;
    (; _cdo, _cds) = cg_con;

    # 8. update fs and fo
    # This is basically equivalent to Kb in Bonan, eq. 14.21
    # TOD reduce allocations
    can_opt.fs      .= _cds ./ cos_sza;
    can_opt.fo      .= _cdo ./ cos_vza;
    can_opt.absfs   .= abs.( can_opt.fs );
    can_opt.absfo   .= abs.( can_opt.fo );
    can_opt.cosΘ_l  .= cos_ttli .* _1s;
    can_opt.cos2Θ_l .= can_opt.cosΘ_l .^ 2;
    can_opt.fsfo    .= can_opt.fs .* can_opt.fo;
    can_opt.absfsfo .= abs.( can_opt.fsfo );

    # 9. probabilities Ps, Po, Pso
    _fac_s = (1 - exp(-ks*Ω*LAI*dx)) / (ks*Ω*LAI*dx);
    _fac_o = (1 - exp(-ko*Ω*LAI*dx)) / (ko*Ω*LAI*dx);
    can_opt.Ps  .= exp.(xl.*ks.*Ω.*LAI) .* _fac_s;
    can_opt.Po  .= exp.(xl.*ko.*Ω.*LAI) .* _fac_o;
    @inline f(x) = psofunction(ko, ks, Ω, LAI, hot, dso, x);

    # TODO minimize the allocations here
    # length(xl) * 7 allocations here!
    @inbounds for j in eachindex(xl)
        can_opt.Pso[j] = quadgk(f, xl[j]-dx, xl[j], rtol=1e-2)[1] / dx;
    end

    # takes care of rounding error
    # can_opt.Pso[can_opt.Pso.>can_opt.Po] =
    # minimum([can_opt.Po[can_opt.Pso.>can_opt.Po]
    #          can_opt.Ps[can_opt.Pso.>can_opt.Po]],dims=2)
    # can_opt.Pso[can_opt.Pso.>can_opt.Ps] =
    # minimum([can_opt.Po[can_opt.Pso.>can_opt.Ps]
    #          can_opt.Ps[can_opt.Pso.>can_opt.Ps]],dims=2)

    return nothing
end

# Canopy matrices

###############################################################################
#
# Update canopy matrices
#
###############################################################################
"""
    canopy_matrices!(leaves::Vector{LeafBios{FT}}, can_opt::CanopyOpticals{FT}) where {FT<:AbstractFloat}

Compute scattering coefficient matrices for direct and diffuse light given
    geometry dependent overall extinction coefficients and pigment dependent
    leaf reflectance and transmission (computed via fluspect). This function
    has to be called before [`short_wave!`](@ref) can be used.
- `leaves` Array of [`LeafBios`](@ref) type struct
- `can_opt` [`CanopyOpticals`](@ref) type struct
"""
function canopy_matrices!(leaves::Vector{LeafBios{FT}}, can_opt::CanopyOpticals{FT}) where {FT<:AbstractFloat}
    # 1. unpack values
    (; ddb, ddf, dob, dof, sdb, sdf, sob, sof) = can_opt;

    # 2. Calculation of reflectance
    nLayer = size(can_opt.sigb)[2];
    @inbounds for i=1:nLayer
        if length(leaves)>1
            τ_SW = leaves[i].τ_SW;
            ρ_SW = leaves[i].ρ_SW;
        else
            τ_SW = leaves[1].τ_SW;
            ρ_SW = leaves[1].ρ_SW;
        end

        #CF: Right now, canopy geometry is the same everywhere,
        #    can be easily extended to layers as well.
        @inbounds for j=1:size(can_opt.sigb, 1)
            # diffuse backscatter coefficient for diffuse incidence
            can_opt.sigb[j,i] = ddb * ρ_SW[j] + ddf * τ_SW[j];
            # diffuse forwardscatter coefficient for diffuse incidence
            can_opt.sigf[j,i] = ddf * ρ_SW[j] + ddb * τ_SW[j];
            # diffuse backscatter coefficient for specular incidence
            can_opt.sb[j,i]   = sdb * ρ_SW[j] + sdf * τ_SW[j];
            # diffuse forwardscatter coefficient for specular incidence
            can_opt.sf[j,i]   = sdf * ρ_SW[j] + sdb * τ_SW[j];
            # directional backscatter  coefficient for diffuse incidence
            can_opt.vb[j,i]   = dob * ρ_SW[j] + dof * τ_SW[j];
            # directional forwardscatter coefficient for diffuse incidence
            can_opt.vf[j,i]   = dof * ρ_SW[j] + dob * τ_SW[j];
            # bidirectional scattering coefficent (directional-directional)
            can_opt.w[j,i]    = sob * ρ_SW[j] + sof * τ_SW[j];
        end
    end

    # 3. attenuation
    can_opt.a .= 1 .- can_opt.sigf;

    return nothing
end

# shortwave

###############################################################################
#
# Simulate short wave radiation
#
###############################################################################
"""
    short_wave!(can::Canopy4RT{FT}, can_opt::CanopyOpticals{FT}, can_rad::CanopyRads{FT}, in_rad::IncomingRadiation{FT}, soil::SoilOpticals{FT}, rt_con::RTCache{FT}) where {FT<:AbstractFloat}

Simulate the short wave radiation through the canopy, given
- `can` [`Canopy4RT`](@ref) type struct
- `can_opt` [`CanopyOpticals`](@ref) type struct
- `can_rad` [`CanopyRads`](@ref) type struct
- `in_rad` [`IncomingRadiation`](@ref) type struct
- `soil` [`SoilOpticals`](@ref) type struct
- `rt_con` [`RTCache`](@ref) type cache
"""
function short_wave!(can::Canopy4RT{FT}, can_opt::CanopyOpticals{FT}, can_rad::CanopyRads{FT}, in_rad::IncomingRadiation{FT}, soil::SoilOpticals{FT}, rt_con::RTCache{FT}) where {FT<:AbstractFloat}
    # unpack values from can and soil
    (; LAI, nLayer, Ω) = can;
    (; ks, sb, sf, sigb) = can_opt;
    (; ρ_SW) = soil;
    sw_con = rt_con.sw_con;

    # 1. define some useful parameters
    iLAI = LAI * Ω / nLayer;

    # 2. scattering and extinction coefficients to
    #    thin layer reflectances and transmittances
    # Eq. 17 in mSCOPE paper (changed here to compute real transmission)
    # this is the original equation: τ_ss = 1 - ks * iLAI;
    τ_ss = exp(-ks * iLAI);
    sw_con.τ_dd .= 1 .- can_opt.a .* iLAI;
    sw_con.τ_sd .= sf   .* iLAI;
    sw_con.ρ_dd .= sigb .* iLAI;
    sw_con.ρ_sd .= sb   .* iLAI;
    (; ρ_dd, ρ_sd, τ_dd, τ_sd) = sw_con;

    # 3. reflectance calculation
    # 3.1 Eq. 18 in mSCOPE paper
    Xss = τ_ss;

    # 3.2 Soil reflectance boundary condition (same for diffuse and direct)
    can_opt.R_sd[:,end] .= ρ_SW;
    can_opt.R_dd[:,end] .= ρ_SW;

    # 3.3 reflectance for each layer from bottom to top
    @inbounds for j in nLayer:-1:1
        sw_con.dnorm      .= 1 .- view(ρ_dd, :, j) .* view(can_opt.R_dd, :, j+1);
        can_opt.Xsd[:,j]  .= ( view(τ_sd, :, j) .+ Xss .* view(can_opt.R_sd, :, j+1) .* view(ρ_dd, :, j) ) ./ sw_con.dnorm;
        can_opt.Xdd[:,j]  .= view(τ_dd, :, j) ./ sw_con.dnorm;
        can_opt.R_sd[:,j] .= view(ρ_sd, :, j) .+ view(τ_dd, :, j) .* ( view(can_opt.R_sd, :, j+1) .* Xss .+ view(can_opt.R_dd, :, j+1) .* view(can_opt.Xsd , :, j) );
        can_opt.R_dd[:,j] .= view(ρ_dd, :, j) .+ view(τ_dd, :, j) .* view(can_opt.R_dd, :, j+1) .* view(can_opt.Xdd , :, j);
    end

    # 4. flux profile calculation
    # Eq. 19 in mSCOPE paper
    # 4.1 Boundary condition at top: Incoming solar radiation
    can_opt.Es_[:,1]    .= in_rad.E_direct;
    can_rad.E_down[:,1] .= in_rad.E_diffuse;

    # 4.2 from top to bottom
    @inbounds for j=1:nLayer
        can_rad.netSW_sunlit[:,j] .= view(can_opt.Es_ , :, j) .* ( 1 .- (τ_ss .+ view(τ_sd, :, j) .+ view(ρ_sd, :, j)) );
        can_opt.Es_[:,j+1]        .= Xss .* view(can_opt.Es_, :, j);
        can_rad.E_down[:,j+1]     .= view(can_opt.Xsd, :, j) .* view(can_opt.Es_, :, j) .+ view(can_opt.Xdd, :, j) .* view(can_rad.E_down, :, j);
        can_rad.E_up[:,j]         .= view(can_opt.R_sd, :, j) .* view(can_opt.Es_, :, j) .+ view(can_opt.R_dd, :, j) .* view(can_rad.E_down, :, j);
    end

    # 4.3 Boundary condition at the bottom, soil reflectance (Lambertian here)
    last_ind_co = lastindex(can_opt.R_sd, 2);
    can_rad.E_up[:,end] .= view(can_opt.R_sd, :, last_ind_co) .* view(can_opt.Es_, :, last_ind_co) .+ view(can_opt.R_dd, :, last_ind_co) .* view(can_rad.E_down, :, last_ind_co);

    # 4.4 Hemispheric total outgoing
    can_rad.Eout .= view(can_rad.E_up, :, 1);

    # 4.5 compute net diffuse radiation per layer:
    @inbounds for j in 1:nLayer
        can_rad.netSW_shade[:,j] .= ( view(can_rad.E_down, :, j) .+ view(can_rad.E_up, :, j+1) ) .* ( 1 .- ( view(τ_dd, :, j) .+ view(ρ_dd, :, j) ) );
        # Add diffuse radiation to direct radiation as well:
        #can_rad.netSW_sunlit[:,j] += can_rad.netSW_shade[:,j]
    end

    # 4.6 outgoing in viewing direction
    # From Canopy
    sw_con.piLoc2 .= can_opt.vb .* view(can_opt.Po, 1:nLayer)' .* view(can_rad.E_down, :, 1:nLayer)  .+
                     can_opt.vf .* view(can_opt.Po, 1:nLayer)' .* view(can_rad.E_up, :, 1:nLayer)  .+
                     can_opt.w  .* view(can_opt.Pso, 1:nLayer)' .* in_rad.E_direct;
    #sw_con.piLoc  .= iLAI .* view(sum(sw_con.piLoc2, dims=2), :, 1);
    @inbounds for j in eachindex(sw_con.piLoc)
        sw_con.piLoc[j] = iLAI * sum( view(sw_con.piLoc2, j, :) );
    end

    # 4.7 From Soil
    sw_con.piLos .= view(can_rad.E_up, :, last_ind_co) .* can_opt.Po[end];
    sw_con.piLo  .= sw_con.piLoc .+ sw_con.piLos;
    can_rad.Lo   .= sw_con.piLo ./ pi;

    # 4.8 Save albedos (hemispheric direct and diffuse and directional (obs))
    # rso and rdo are not computed separately
    can_rad.alb_obs     .= sw_con.piLo ./ ( in_rad.E_direct .+ in_rad.E_diffuse );
    can_rad.alb_direct  .= view(can_opt.R_sd, :, 1);
    can_rad.alb_diffuse .= view(can_opt.R_dd, :, 1);

    return nothing
end

# Canopy fluxes

"""
    canopy_fluxes!(
                can::Canopy4RT{FT},
                can_opt::CanopyOpticals{FT},
                can_rad::CanopyRads{FT},
                in_rad::IncomingRadiation{FT},
                soil::SoilOpticals{FT},
                leaves::Vector{LeafBios{FT}},
                wls::WaveLengths{FT},
                rt_con::RTCache{FT}
    ) where {FT<:AbstractFloat}

Computes a variety of integrated fluxes from the spectrally resolved
    computations in the short-wave Canopy RT (e.g. absorbed soil radiation,
    absorbed direct and diffuse PAR by layer (and angles for direct), net
    direct and diffuse energy balance per layer), given
- `can` [`Canopy4RT`](@ref) type struct
- `can_opt` [`CanopyOpticals`](@ref) type struct
- `can_rad` [`CanopyRads`](@ref) type struct
- `in_rad` [`IncomingRadiation`](@ref) type struct
- `soil` [`SoilOpticals`](@ref) type struct
- `leaves` Array of [`LeafBios`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
- `rt_con` [`RTCache`](@ref) type cache
"""
function canopy_fluxes!(
            can::Canopy4RT{FT},
            can_opt::CanopyOpticals{FT},
            can_rad::CanopyRads{FT},
            in_rad::IncomingRadiation{FT},
            soil::SoilOpticals{FT},
            leaves::Vector{LeafBios{FT}},
            wls::WaveLengths{FT},
            rt_con::RTCache{FT}
) where {FT<:AbstractFloat}
    # 1. unpack variables from structures
    (; LAI, nLayer) = can;
    (; ε_SW) = soil;
    (; dWL, dWL_iPAR, dWL_iPAR_700, iPAR, iPAR_700, WL_iPAR) = wls;
    cf_con = rt_con.cf_con;

    # 2. compute some useful variables
    tLAI = LAI / nLayer;
    fac  = FT(1e-3);

    # 3. Compute some fluxes, can be done separately if needed
    #    this is absolute fluxes now, for the entire soil
    last_ind_cr            = lastindex(can_rad.E_down,2);
    cf_con.abs_wave       .= view(can_rad.E_down, :, last_ind_cr) .* ε_SW;
    can_rad.RnSoil_diffuse = fac * numerical∫(cf_con.abs_wave, dWL);
    cf_con.abs_wave       .= view(can_opt.Es_, :, last_ind_cr) .* ε_SW;
    can_rad.RnSoil_direct  = fac * numerical∫(cf_con.abs_wave, dWL);
    can_rad.RnSoil         = can_rad.RnSoil_direct + can_rad.RnSoil_diffuse;

    # 4. Normalization factor for leaf direct PAR
    #    weighted sum has to be 1 to conserve net SW direct
    #    Direct PAR is normalized by layer Ps value
    mul!(cf_con.absfs_lidf, adjoint(can_opt.absfs), can.lidf);
    normi       = 1 / mean(cf_con.absfs_lidf);
    cf_con.lPs .= (view(can_opt.Ps, 1:nLayer) .+ view(can_opt.Ps, 2:nLayer+1)) ./ 2;
    (; lPs) = cf_con;

    @inbounds for j in 1:nLayer
        if length(leaves)>1
            cf_con.kChlrel .= view(leaves[j].kChlrel, iPAR);
        else
            cf_con.kChlrel .= view(leaves[1].kChlrel, iPAR);
        end

        # for diffuse PAR
        cf_con.E_iPAR .= view(can_rad.netSW_shade, iPAR, j);
        e2phot!(WL_iPAR, cf_con.E_iPAR, cf_con.PAR_diff);
        cf_con.PAR_diff .*= fac / tLAI;

        # for direct PAR
        cf_con.E_iPAR .= view(can_rad.netSW_sunlit, iPAR, j);
        e2phot!(WL_iPAR, cf_con.E_iPAR, cf_con.PAR_dir);
        cf_con.PAR_dir .*= fac / tLAI / lPs[j];

        # for leaf absorbed
        cf_con.PAR_diffCab .= cf_con.kChlrel .* cf_con.PAR_diff;
        cf_con.PAR_dirCab  .= cf_con.kChlrel .* cf_con.PAR_dir;

        # Absorbed PAR per leaf for shaded, PAR_diff and PAR_dir changed!
        _dif = numerical∫(cf_con.PAR_diff, dWL_iPAR);
        _dir = numerical∫(cf_con.PAR_dir, dWL_iPAR) * normi;
        can_rad.absPAR_shade[j] = _dif;
        can_rad.absPAR_sun[:,:,j] .= can_opt.absfs .* _dir;

        # absorbed PAR for photosynthesis (set it to be the lower of 2*X_700 or X_750)
        _difCab_700 = numerical∫(cf_con.PAR_diffCab[iPAR_700], dWL_iPAR_700);
        _dirCab_700 = numerical∫(cf_con.PAR_dirCab[iPAR_700] , dWL_iPAR_700) * normi;
        _difCab_750 = numerical∫(cf_con.PAR_diffCab, dWL_iPAR);
        _dirCab_750 = numerical∫(cf_con.PAR_dirCab , dWL_iPAR) * normi;
        _difCab = min(2*_difCab_700, _difCab_750);
        _dirCab = min(2*_dirCab_700, _dirCab_750);
        can_rad.absPAR_shadeCab[j] = _difCab;
        can_rad.absPAR_sunCab[:,:,j] .= can_opt.absfs .* _dirCab;
        can_rad.absPAR_sunCab[:,:,j] .+= _difCab;
    end

    # 5. Total PAR
    # TODO considering remove this part, if we are not using it
    #    re-use the PAR_dir and PAR_diff in the rt_con
    cf_con.E_iPAR .= view(in_rad.E_direct, iPAR);
    e2phot!(WL_iPAR, cf_con.E_iPAR, cf_con.PAR_dir);
    cf_con.E_iPAR .= view(in_rad.E_diffuse, iPAR);
    e2phot!(WL_iPAR, cf_con.E_iPAR, cf_con.PAR_diff);
    can_rad.incomingPAR_direct = fac * numerical∫(cf_con.PAR_dir , dWL_iPAR);
    can_rad.incomingPAR_diffuse = fac * numerical∫(cf_con.PAR_diff, dWL_iPAR);
    can_rad.incomingPAR = can_rad.incomingPAR_diffuse + can_rad.incomingPAR_direct;
    @inbounds for i in 1:nLayer
        cf_con.E_all .= view(can_rad.netSW_shade, :, i);
        can_rad.intNetSW_shade[i] = numerical∫(cf_con.E_all, dWL) * fac / tLAI;
        cf_con.E_all .= view(can_rad.netSW_sunlit, :, i);
        can_rad.intNetSW_sunlit[i] = numerical∫(cf_con.E_all, dWL) * fac / tLAI / lPs[i] + can_rad.intNetSW_shade[i];
    end

    return nothing
end

# SIF flux

###############################################################################
#
# Calculate SIF fluxes
#
###############################################################################
"""
    SIF_fluxes!(leaves::Vector{LeafBios{FT}},
                can_opt::CanopyOpticals{FT},
                can_rad::CanopyRads{FT},
                can::Canopy4RT{FT},
                soil::SoilOpticals{FT},
                wls::WaveLengths{FT},
                rt_con::RTCache{FT},
                rt_dim::RTDimensions;
                photon::Bool = true
    ) where {FT<:AbstractFloat}

Computes 2-stream diffusive radiation transport for SIF radiation (calls
    [`diffusive_S!`] internally). Layer reflectance and transmission is
    computed from SW optical properties, layer sources from absorbed light and
    SIF efficiencies. Boundary conditions are zero SIF incoming from atmosphere
    or soil.
- `leaves` Array of [`LeafBios`](@ref) type struct
- `can_opt` [`CanopyOpticals`](@ref) type struct
- `can_rad` [`CanopyRads`](@ref) type struct
- `can` [`Canopy4RT`](@ref) type struct
- `soil` [`SoilOpticals`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
- `rt_con` [`RTCache`](@ref) type cache
- `rt_dim` [`RTDimensions`](@ref) type struct
- `photon` If true, use photon unit in the matrix conversion

"""
function SIF_fluxes!(
            leaves::Vector{LeafBios{FT}},
            can_opt::CanopyOpticals{FT},
            can_rad::CanopyRads{FT},
            can::Canopy4RT{FT},
            soil::SoilOpticals{FT},
            wls::WaveLengths{FT},
            rt_con::RTCache{FT},
            rt_dim::RTDimensions;
            photon::Bool = true
) where {FT<:AbstractFloat}
    # unpack variables from structures
    (; LAI, lidf, nLayer, Ω) = can;
    (; a, absfo, absfs, absfsfo, cosΘ_l, cos2Θ_l, fo, fs, fsfo, Po, Ps, Pso, sigb, vb, vf) = can_opt;
    (; E_down, E_up, ϕ_shade, ϕ_sun) = can_rad;
    (; ρ_SW_SIF) = soil;
    (; dWL_iWLE, iWLE, iWLF, WLE, WLF) = wls;
    sf_con = rt_con.sf_con;

    # 1. define some useful parameters
    iLAI = LAI * Ω / nLayer;

    # 2. calculate some useful parameters
    sf_con.τ_dd .= 1 .- view(a, iWLF, :) .* iLAI;
    sf_con.ρ_dd .= view(sigb, iWLF, :) .* iLAI;

    # 3. Compute mid layer Ps,Po,Pso
    #    Qso = (Pso[1:nLayer] + Pso[2:nLayer+1]) / 2;
    #    Qs  = ( Ps[1:nLayer] +  Ps[2:nLayer+1]) / 2;
    #    Qo  = ( Po[1:nLayer] +  Po[2:nLayer+1]) / 2;
    Qso = view(Pso, 1:nLayer);
    Qs  = view(Ps , 1:nLayer);
    Qo  = view(Po , 1:nLayer);

    # 4.  reflectance, transmittance factors in a thin layer the following are
    #     vectors with length [nl,nWL]
    sf_con.sun_dwl_iWlE .= view(can_opt.Es_, iWLE, 1) .* dWL_iWLE;
    if photon sf_con.sun_dwl_iWlE .*= WLE .* _FAC(FT) end;
    @inbounds for i=1:nLayer
        if length(leaves)>1
            Mb = leaves[i].Mb;
            Mf = leaves[i].Mf;
        else
            Mb = leaves[1].Mb;
            Mf = leaves[1].Mf;
        end
        sf_con.M⁺ .= (Mb .+ Mf) ./ 2;
        sf_con.M⁻ .= (Mb .- Mf) ./ 2;

        # Need to normalize incoming radiation bin
        # change from mSCOPE to enable different WL grids!
        mul!(sf_con.M⁻_sun, sf_con.M⁻, sf_con.sun_dwl_iWlE);
        mul!(sf_con.M⁺_sun, sf_con.M⁺, sf_con.sun_dwl_iWlE);

        sf_con.tmp_dwl_iWlE .= view(E_down, iWLE, i) .* dWL_iWLE;
        if photon sf_con.tmp_dwl_iWlE .*= WLE .* _FAC(FT) end;
        mul!(sf_con.M⁺⁻, sf_con.M⁺, sf_con.tmp_dwl_iWlE);

        sf_con.tmp_dwl_iWlE .= view(E_up, iWLE, i+1) .* dWL_iWLE;
        if photon sf_con.tmp_dwl_iWlE .*= WLE .* _FAC(FT) end;
        mul!(sf_con.M⁺⁺, sf_con.M⁺, sf_con.tmp_dwl_iWlE);

        sf_con.tmp_dwl_iWlE .= view(E_up, iWLE, i+1) .* dWL_iWLE;
        if photon sf_con.tmp_dwl_iWlE .*= WLE .* _FAC(FT) end;
        mul!(sf_con.M⁻⁺, sf_con.M⁻, sf_con.tmp_dwl_iWlE);

        sf_con.tmp_dwl_iWlE .= view(E_down, iWLE, i) .* dWL_iWLE;
        if photon sf_con.tmp_dwl_iWlE .*= WLE .* _FAC(FT) end;
        mul!(sf_con.M⁻⁻, sf_con.M⁻, sf_con.tmp_dwl_iWlE);

        # Here comes the tedious part:
        # TODO move them to a seprate function
        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* cosΘ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        sunCos = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= ϕ_shade[i] .* cosΘ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        shadeCos = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* cos2Θ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        sunCos2 = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= ϕ_shade[i] .* cos2Θ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        shadeCos2 = mean(sf_con.ϕ_cosΘ_lidf);

        mul!(sf_con.ϕ_cosΘ_lidf, view(ϕ_sun, :, :, i)', lidf);
        sunLidf = mean(sf_con.ϕ_cosΘ_lidf);
        shadeLidf = mean(lidf) * ϕ_shade[i];

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* absfsfo;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_absfsfo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* fsfo;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_fsfo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.wfEs .= _mean_absfsfo .* sf_con.M⁺_sun .+ _mean_fsfo .* sf_con.M⁻_sun;

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* absfs;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_absfs = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* fs .* cosΘ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_fs = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.sfEs .= _mean_absfs .* sf_con.M⁺_sun .- _mean_fs .* sf_con.M⁻_sun;
        sf_con.sbEs .= _mean_absfs .* sf_con.M⁺_sun .+ _mean_fs .* sf_con.M⁻_sun;

        sf_con.ϕ_cosΘ .= ϕ_shade[i] .* absfo;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_absfo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= ϕ_shade[i] .* fo .* cosΘ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_fo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.vfEplu_shade .= _mean_absfo .* sf_con.M⁺⁺ .- _mean_fo .* sf_con.M⁻⁺;
        sf_con.vbEmin_shade .= _mean_absfo .* sf_con.M⁺⁻ .+ _mean_fo .* sf_con.M⁻⁻;

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* absfo;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_absfo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* fo .* cosΘ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_fo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.vfEplu_sun .= _mean_absfo .* sf_con.M⁺⁺ .- _mean_fo .* sf_con.M⁻⁺;
        sf_con.vbEmin_sun .= _mean_absfo .* sf_con.M⁺⁻ .+ _mean_fo .* sf_con.M⁻⁻;

        sf_con.sigfEmin_shade .= shadeLidf .* sf_con.M⁺⁻ .- shadeCos2 .* sf_con.M⁻⁻;
        sf_con.sigbEmin_shade .= shadeLidf .* sf_con.M⁺⁻ .+ shadeCos2 .* sf_con.M⁻⁻;
        sf_con.sigfEmin_sun   .= sunLidf   .* sf_con.M⁺⁻ .- sunCos2   .* sf_con.M⁻⁻;
        sf_con.sigbEmin_sun   .= sunLidf   .* sf_con.M⁺⁻ .+ sunCos2   .* sf_con.M⁻⁻;
        sf_con.sigfEplu_shade .= shadeLidf .* sf_con.M⁺⁺ .- shadeCos2 .* sf_con.M⁻⁺;
        sf_con.sigbEplu_shade .= shadeLidf .* sf_con.M⁺⁺ .+ shadeCos2 .* sf_con.M⁻⁺;
        sf_con.sigfEplu_sun   .= sunLidf   .* sf_con.M⁺⁺ .- sunCos2   .* sf_con.M⁻⁺;
        sf_con.sigbEplu_sun   .= sunLidf   .* sf_con.M⁺⁺ .+ sunCos2   .* sf_con.M⁻⁺;

        # Fluxes:
        #    sunlit for each layer
        #    shade leaf for each layer
        #    Eq. 29a for sunlit leaf
        #    Eq. 29b for sunlit leaf
        #    Eq. 29a for shade leaf
        #    Eq. 29b for shade leaf
        sf_con.piLs[:,i]  .= sf_con.wfEs .+ sf_con.vfEplu_sun .+ sf_con.vbEmin_sun;
        sf_con.piLd[:,i]  .= sf_con.vbEmin_shade .+ sf_con.vfEplu_shade;
        sf_con.Fsmin[:,i] .= sf_con.sfEs .+ sf_con.sigfEmin_sun .+ sf_con.sigbEplu_sun;
        sf_con.Fsplu[:,i] .= sf_con.sbEs .+ sf_con.sigbEmin_sun .+ sf_con.sigfEplu_sun;
        sf_con.Fdmin[:,i] .= sf_con.sigfEmin_shade .+ sf_con.sigbEplu_shade;
        sf_con.Fdplu[:,i] .= sf_con.sigbEmin_shade .+ sf_con.sigfEplu_shade;

        # Total weighted fluxes
        _qs_iLAI   = Qs[i] * iLAI;
        _1_qs_iLAI = (1 - Qs[i]) * iLAI;
        sf_con.S⁻[:,i]   .= _qs_iLAI .* view(sf_con.Fsmin, :, i) .+ _1_qs_iLAI .* view(sf_con.Fdmin, :, i);
        sf_con.S⁺[:,i]   .= _qs_iLAI .* view(sf_con.Fsplu, :, i) .+ _1_qs_iLAI .* view(sf_con.Fdplu, :, i);
        sf_con.Femo[:,i] .= _qs_iLAI .* view(sf_con.piLs , :, i) .+ _1_qs_iLAI .* view(sf_con.piLd , :, i);
    end

    # 5. Compute diffusive fluxes within canopy
    #    Use Zero SIF fluxes as top and bottom boundary:
    # TODO pre-allocate these!
    diffusive_S!(sf_con, soil, rt_dim);

    # 6. Save in output structures!
    #    direct Sunlit leaves
    #    direct shaded leaves
    _iLAI_pi = iLAI / FT(pi);

    # why this step is so slow?
    sf_con.tmp_1d_nLayer .= Qso .* _iLAI_pi;
    mul!(can_rad.SIF_obs_sunlit, sf_con.piLs, sf_con.tmp_1d_nLayer);
    sf_con.tmp_1d_nLayer .= (Qo .- Qso) .* _iLAI_pi;
    mul!(can_rad.SIF_obs_shaded, sf_con.piLd, sf_con.tmp_1d_nLayer);

    # 7. SIF from scattered internally and soil contribution
    sf_con.tmp_2d_nWlF_nLayer .= view(vb, iWLF, :) .* view(sf_con.F⁻, :, 1:nLayer) .+ view(vf, iWLF, :) .* view(sf_con.F⁺, :, 1:nLayer);
    mul!(can_rad.SIF_obs_scattered, sf_con.tmp_2d_nWlF_nLayer, Qo);
    can_rad.SIF_obs_scattered .*= _iLAI_pi;
    can_rad.SIF_obs_soil .= ( ρ_SW_SIF .* view(sf_con.F⁻, :, nLayer+1) ) .* Po[end] ./ FT(pi);

    can_rad.SIF_hemi .= view(sf_con.F⁺, :, 1);
    can_rad.SIF_obs  .= can_rad.SIF_obs_sunlit .+ can_rad.SIF_obs_shaded .+ can_rad.SIF_obs_scattered .+ can_rad.SIF_obs_soil;
    # can_rad.SIF_sum[:]  = sum(sf_con.S⁻ + sf_con.S⁺, dims=2);
    @inbounds for j in eachindex(can_rad.SIF_sum)
        can_rad.SIF_sum[j] = sum( view(sf_con.S⁻, j, :) ) + sum( view(sf_con.S⁺, j, :) );
    end

    if photon
        can_rad.SIF_obs_sunlit ./= WLF .* _FAC(FT);
        can_rad.SIF_obs_shaded ./= WLF .* _FAC(FT);
        can_rad.SIF_obs_scattered ./= WLF .* _FAC(FT);
        can_rad.SIF_obs_soil ./= WLF .* _FAC(FT);
        can_rad.SIF_obs ./= WLF .* _FAC(FT);
    end;

    return nothing
end


"""
    SIF_fluxes!(leaf::LeafBios{FT}, in_rad::IncomingRadiation{FT}, wls::WaveLengths{FT}, rt_con::RTCache{FT}, fqe::FT = FT(0.01); photon::Bool = true) where {FT<:AbstractFloat}

Leaf level SIF, given
- `leaf` [`LeafBios`](@ref) type struct
- `in_rad` [`IncomingRadiation`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
- `rt_con` [`RTCache`](@ref) type cache
- `fqe` Fluorescence quantum yield (default at 1%)
- `photon` If true, use photon unit in the matrix conversion

Note that `in_rad` assumes direct light with zenith angle of 0, and a zenith
    angle correction needs to be made before passing it to this function. The
    up- and down-ward SIF are stored in `sf_con` as `M⁻_sun` and `M⁺_sun`.
"""
function SIF_fluxes!(leaf::LeafBios{FT}, in_rad::IncomingRadiation{FT}, wls::WaveLengths{FT}, rt_con::RTCache{FT}, fqe::FT = FT(0.01); photon::Bool = true) where {FT<:AbstractFloat}
    # unpack the values
    (; Mb, Mf) = leaf;
    (; dWL_iWLE, iWLE, WLE, WLF) = wls;
    sf_con = rt_con.sf_con;
    sf_con.tmp_dwl_iWlE  .= (view(in_rad.E_direct , iWLE, 1) .+ view(in_rad.E_diffuse, iWLE, 1)) .* dWL_iWLE;
    if photon
        sf_con.tmp_dwl_iWlE .*= WLE .* _FAC(FT);
    end;

    # calculate the SIF spectra for direct light
    # sf_con.M⁺ .= (Mb .+ Mf) ./ 2;
    # sf_con.M⁻ .= (Mb .- Mf) ./ 2;
    # mul!(sf_con.M⁺_sun, sf_con.M⁺, sf_con.tmp_dwl_iWlE);
    # mul!(sf_con.M⁻_sun, sf_con.M⁻, sf_con.tmp_dwl_iWlE);
    mul!(sf_con.M⁺_sun, Mb, sf_con.tmp_dwl_iWlE);
    mul!(sf_con.M⁻_sun, Mf, sf_con.tmp_dwl_iWlE);
    if photon
        sf_con.M⁺_sun ./= WLF .* _FAC(FT);
        sf_con.M⁻_sun ./= WLF .* _FAC(FT);
    end;

    # divide by pi to account for scattering
    sf_con.M⁻_sun .*= fqe / pi;
    sf_con.M⁺_sun .*= fqe / pi;

    return nothing
end