#=
using Revise; include(joinpath("test", "collection", "expr_code_lowered_single_expression.jl"))
=#
using Test
import MultiBroadcastFusion as MBF

@testset "code_lowered_single_expression" begin
    expr_in = :(@. y1 = x1 + x2 + x3 + x4)
    expr_out = :(Base.materialize!(y1, Base.broadcasted(+, x1, x2, x3, x4)))
    @test MBF.code_lowered_single_expression(expr_in) == expr_out
end

@testset "code_lowered_single_expression - examples from the wild" begin
    expr_in = quote
        @. ᶜcloud_fraction = quad_loop(
            SG_quad,
            ᶜts,
            Geometry.WVector(p.precomputed.ᶜgradᵥ_q_tot),
            Geometry.WVector(p.precomputed.ᶜgradᵥ_θ_liq_ice),
            coeff,
            ᶜmixing_length,
            thermo_params,
        )
    end
    expr_out = :(Base.materialize!(
        ᶜcloud_fraction,
        Base.broadcasted(
            quad_loop,
            SG_quad,
            ᶜts,
            Base.broadcasted(
                Base.getproperty(Geometry, :WVector),
                Base.getproperty(
                    Base.getproperty(p, :precomputed),
                    :ᶜgradᵥ_q_tot,
                ),
            ),
            Base.broadcasted(
                Base.getproperty(Geometry, :WVector),
                Base.getproperty(
                    Base.getproperty(p, :precomputed),
                    :ᶜgradᵥ_θ_liq_ice,
                ),
            ),
            coeff,
            ᶜmixing_length,
            thermo_params,
        ),
    ))
    @test MBF.code_lowered_single_expression(expr_in) == expr_out
end
