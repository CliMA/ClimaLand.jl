#=
propertynames(p[:soil])
(:K, :ψ, :θ_l, :T, :κ, :turbulent_fluxes, :R_n, :top_bc, :infiltration, :bottom_bc)
propertynames(p[:soilco2])
(:D, :Sm, :top_bc, :bottom_bc)
propertynames(p[:canopy])
(:hydraulics, :conductance, :photosynthesis, :radiative_transfer, :autotrophic_respiration, :energy)
propertynames(p[:drivers])
(:P_liq, :P_snow, :T, :P, :u, :q, :c_co2, :SW_d, :LW_d, :θs)
propertynames(p[:soil][:turbulent_fluxes])
(:lhf, :shf, :vapor_flux, :r_ae)
propertynames(p[:canopy][:hydraulics])
(:β, :ψ, :fa, :fa_roots, :area_index)
propertynames(p[:canopy][:conductance])
(:medlyn_term, :gs, :transpiration)
propertynames(p[:canopy][:photosynthesis])
(:An, :GPP, :Rd, :Vcmax25)
propertynames(p[:canopy][:radiative_transfer])
(:apar, :par, :rpar, :tpar, :anir, :nir, :rnir, :tnir, :LW_n, :SW_n)
propertynames(p[:canopy][:autotrophic_respiration])
(:Ra,)
propertynames(p[:canopy][:energy])
(:shf, :lhf, :fa_energy_roots, :r_ae)
propertynames(Y.soilco2)
(:C,)
propertynames(Y.soil)
(:ϑ_l, :θ_i, :ρe_int)
propertynames(Y.canopy)
(:hydraulics, :energy)
=#

###
# GPP
###
add_diagnostic_variable!(
    short_name = "GPP",
    long_name = "Gross Primary Productivity",
    standard_name = "GPP",
    units = "gC m^-2 s^-1",
    compute! = (out, Y, p, t) -> begin
        if isnothing(out)
            return copy(p.canopy.photosynthesis.GPP)
        else
            out .= p.canopy.photosynthesis.GPP
        end
    end,
)





















