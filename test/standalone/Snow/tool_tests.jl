using ClimaLand
using Test

import ClimaParams as CP
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using Thermodynamics
using SurfaceFluxes
using StaticArrays
using Dates
using ClimaLand.Snow
using ClimaLand.Domains
import ClimaLand.Parameters as LP
using ClimaComms
ClimaComms.@import_required_backends

using Flux, Adapt, JLD2, InteractiveUtils #leftover weak dep for ConstrainedNeuralModelExt
using Downloads, Statistics, DataFrames #leftover weak dep for SNOTELScraperExt


SNOTELScraperExt = Base.get_extension(ClimaLand, :SNOTELScraperExt)
ConstrainedNeuralModelExt =
    Base.get_extension(ClimaLand, :ConstrainedNeuralModelExt)
NeuralSnowExt = Base.get_extension(ClimaLand, :ConstrainedNeuralModelExt)

if !isnothing(ConstrainedNeuralModelExt)
    CNM = ConstrainedNeuralModelExt.ConstrainedNeuralModels
    import .ConstrainedNeuralModelExt.ConstrainedNeuralModels.@bound as @bound
    import .ConstrainedNeuralModelExt.ConstrainedNeuralModels.@bound_type as @bound_type
    mt = SMatrix{1, N, <:AbstractFloat, N} where {N <: Int} #to deal with one eval() call in testing
end

if !isnothing(SNOTELScraperExt)
    DataTools = SNOTELScraperExt.DataTools
    NeuralSnow = ConstrainedNeuralModelExt.NeuralSnow
    @testset "Testing SNOTELScraper Extension" begin
        start_date = "2015-01-01"
        end_date = "2015-02-01"
        station_id = 1030
        station_state = "CO"
        real_link = "https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customSingleStationReport/daily/start_of_period/1030:CO:SNTL/2015-01-01,2015-02-01/SNWD::value,"
        @test DataTools.data_url(
            station_id,
            station_state,
            ["SNWD"],
            start_date = start_date,
            end_date = end_date,
        ) == real_link

        colnames = [
            "date",
            "SWE",
            "z",
            "precip",
            "rel_hum_avg",
            "sol_rad_avg",
            "wind_speed_avg",
            "air_temp_avg",
        ]
        # Leave out tests that scrape data, to limit failed builds due to SNOTEL site being down:
        #=
        test_data1 = DataTools.get_data(
            station_id,
            station_state,
            ["SNWD"],
            start = start_date,
            finish = end_date,
        )
        @test typeof(test_data1) == DataFrame
        @test size(test_data1) == (32, 2)
        @test typeof(test_data1[1, 1]) == Date
        @test sum(test_data1[!, 2]) == 1168
        =#

        download_link = "https://caltech.box.com/shared/static/4tih9hiydrc7bcvrpkr2x727b5l54oiq.csv"
        download_data = DataTools.df_from_url(download_link) #CSV.read(HTTP.get(download_link).body, DataFrame)
        test_data2 = deepcopy(download_data)
        #=
        test_data2 = DataTools.sitedata_daily(
            station_id,
            station_state,
            start = start_date,
            finish = end_date,
        )
        =#
        @test typeof(test_data2) === DataFrame
        @test size(test_data2) == (32, 8)
        @test DataFrames.names(test_data2) == colnames
        @test typeof(test_data2[1, 1]) == Date
        @test sum(test_data2[!, :z]) == 1168
        #@test isequal(download_data, test_data2)

        download_link_hr = "https://caltech.box.com/shared/static/4q7qtnc2nv8rwh4q8s0afe238ydzfl93.csv"
        download_data = DataTools.df_from_url(download_link_hr) #CSV.read(HTTP.get(download_link_hr).body, DataFrame)
        test_data3 = deepcopy(download_data)
        #=
        test_data3 = DataTools.sitedata_hourly(
            station_id,
            station_state,
            start = start_date,
            finish = end_date,
        )
        =#
        test_data3[!, :id] .= station_id
        test_data2[!, :id] .= station_id
        push!(colnames, "id")
        @test typeof(test_data3) === DataFrame
        @test size(test_data3) == (768, 9)
        @test DataFrames.names(test_data3) == colnames
        @test typeof(test_data3[1, 1]) == DateTime
        @test sum(skipmissing(test_data3[!, :z])) == 25063
        @test sum(skipmissing(test_data3[!, :rel_hum_avg])) == 62232
        #@test isequal(download_data, test_data3)

        test_bounds = Dict{Symbol, Tuple{Real, Real}}(
            :z => (0, 37),
            :rel_hum_avg => (84, 100),
        )
        test_data4 = DataTools.apply_bounds(test_data3, test_bounds)
        @test sum(skipmissing(test_data4[!, :z])) == 19007
        @test sum(skipmissing(test_data4[!, :rel_hum_avg])) == 45973

        test_data5 = DataTools.hourly2daily(test_data3)
        @test sum(test_data5[!, :z]) == 1168
        @test sum(test_data5[!, :rel_hum_avg]) == 2593
        @test sum(describe(test_data5, :nmissing)[!, 2]) == 0

        test_data6 = DataTools.rectify_daily_hourly(test_data2, test_data5)
        @test sum(test_data6[!, :z]) == 1168
        @test sum(test_data6[!, :rel_hum_avg]) == 2593
        inch2meter = 0.0254
        kmphr2mps = 5.0 / 18.0
        scales = Dict{Symbol, Real}(
            :SWE => inch2meter,
            :z => inch2meter,
            :precip => inch2meter,
            :rel_hum_avg => 0.01,
            :wind_speed_avg => kmphr2mps,
        )
        test_data7 = DataTools.scale_cols(test_data6, scales)
        @test sum(test_data7[!, :rel_hum_avg]) ≈ 25.93 atol = 1e-3
        test_data8 = DataTools.makediffs(test_data7, Day(1))
        @test size(test_data8) == (31, 12)
        @test mean(test_data8[!, :dzdt]) ≈ 0 atol = 1e-7
        @test mean(test_data8[!, :dSWEdt]) ≈ 0 atol = 1e-7
        @test minimum(test_data8[!, :dprecipdt]) == 0

        test_data9 = DataTools.rolldata(test_data8, Day(1), 7)
        @test size(test_data9) == (25, 12)
        @test mean(test_data9[!, :dzdt]) ≈ 0 atol = 1e-8
        @test mean(test_data9[!, :dSWEdt]) ≈ 0 atol = 1e-8
        @test minimum(test_data9[!, :dprecipdt]) == 0
        @test minimum(test_data9[!, :z]) ≈ 0.8636 atol = 1e-4

        test_data10 = DataTools.prep_data(
            test_data8,
            extract_vars = [:dprecipdt_rain, :dprecipdt_snow],
        )
        x_train, y_train = DataTools.make_data(
            test_data10,
            [:dprecipdt_rain],
            :dprecipdt_snow,
            1.0,
        )
        @test test_data10[!, 1] .+ test_data10[!, 2] ==
              test_data8[!, :dprecipdt]
        @test size(x_train) == (1, 31)
        @test size(y_train) == (1, 31)

        @test DataTools.fix_temp([10.0], alaska = true)[1] ≈ 9.968 atol = 1e-3

        x = ones(366, 7)
        x = Float32.(x)
        x = DataFrame(
            x,
            [
                :SWE,
                :z,
                :precip,
                :rel_hum_avg,
                :sol_rad_avg,
                :wind_speed_avg,
                :air_temp_avg,
            ],
        )
        x[!, :id] .= station_id
        x[!, :date] .= collect(Date("2015-10-01"):Day(1):Date("2016-09-30"))
        x[200, :precip] = 2.0
        x[49, :precip] = 2.0
        x[49, :SWE] = 3.0
        allowmissing!(x)
        x[[3, 4, 10], :sol_rad_avg] .= missing
        x[100, :wind_speed_avg] = 1000.0
        #x1 = DataTools.serreze_qc(x, station_id, station_state)
        #@test size(x1) == (365, 9)
        #@test count(ismissing, x1[!, :air_temp_avg]) == 365
        #@test sum(x1[!, :precip]) == 365
        x1, n1 = DataTools.yan_qc(x)
        @test count(ismissing, x1[!, :SWE]) == 366
        @test n1 == 1.5
        x1 = DataTools.livneh_bc(x)
        @test sum(x1[!, :precip]) == 852
        @test sum(x1[!, :SWE]) == 368
        x1 = DataTools.z_qc(x)
        @test sum(x1) == 153
        x1, x2, x3 = DataTools.manual_filter(x)
        @test sum(x1) == 0
        @test sum(x2) == 30
        @test sum(x3) == 0
        x1 = DataTools.impute_data(x, t1 = 1, t2 = 2, dt = Day(1))
        @test sum(x1) == 366
        x1 = DataTools.qc_filter(x)
        @test sum(x1) == 1
    end

    @testset "Testing ConstrainedNeuralModels module" begin
        FT = Float32
        no_scale = CNM.NoScaling{FT}()
        use_in_scales = ones(FT, 7)
        c_scale = CNM.ConstScaling(use_in_scales)
        @test typeof(no_scale) <: CNM.InputWeighting{FT}
        @test typeof(c_scale) <: CNM.InputWeighting{FT}
        @test c_scale.in_scales == use_in_scales

        testfunc1(x, y) = x .^ 2 .+ y
        ub = CNM.UpperOnly{typeof(testfunc1)}(testfunc1)
        lb = CNM.LowerOnly{typeof(testfunc1)}(testfunc1)
        @test typeof(ub) <: CNM.ConstraintType
        @test typeof(lb) <: CNM.ConstraintType
        @test CNM.output_dim(ub) == 2
        @test CNM.output_dim(lb) == 2
        @test ub.bound(3, 1) == 10
        @test lb.bound(3, 1) == 10
        @test ub.bound == testfunc1
        @test lb.bound == testfunc1
        ts = CNM.TwoSided{typeof(testfunc1), typeof(testfunc1)}(
            testfunc1,
            testfunc1,
        )
        @test typeof(ts) <: CNM.ConstraintType
        @test CNM.output_dim(ts) == 3
        @test ts.upper_bound(3, 1) == 10
        @test ts.lower_bound(3, 1) == 10
        @test ts.upper_bound == testfunc1
        @test ts.lower_bound == testfunc1
        @test CNM.get_bounds(ub, 4, 3) == (19,)
        @test CNM.get_bounds(lb, 4, 3) == (19,)
        @test CNM.get_bounds(ts, 4, 3) == (19, 19)
        @test CNM.get_bounds(ub) == (testfunc1,)
        @test CNM.get_bounds(lb) == (testfunc1,)
        @test CNM.get_bounds(ts) == (testfunc1, testfunc1)

        inp_1 = Matrix{FT}([1.0 2.0])
        ans_1 = [2 6; 1 2]
        ans_11 = [2 6; 2 6; 1 2]
        ub_1 = CNM.boundary_connection(ub, inp_1, inp_1)
        lb_1 = CNM.boundary_connection(lb, inp_1, inp_1)
        ts_1 = CNM.boundary_connection(ts, inp_1, inp_1)

        @test ub_1 == ans_1
        @test lb_1 == ans_1
        @test ts_1 == ans_11
        @test typeof(ub_1) == Matrix{FT}
        @test typeof(lb_1) == Matrix{FT}
        @test typeof(ts_1) == Matrix{FT}
        @test eltype(ub_1) == FT
        @test eltype(lb_1) == FT
        @test eltype(ts_1) == FT

        inp_2 = SMatrix{1, 2, FT}(inp_1)
        ans_2 = SMatrix{2, 2, FT}(ans_1)
        ans_22 = SMatrix{3, 2, FT}(ans_11)
        ub_2 = CNM.boundary_connection(ub, inp_2, inp_2)
        lb_2 = CNM.boundary_connection(lb, inp_2, inp_2)
        ts_2 = CNM.boundary_connection(ts, inp_2, inp_2)

        @test ub_2 == ans_2
        @test lb_2 == ans_2
        @test ts_2 == ans_22
        @test typeof(ub_2) == SMatrix{2, 2, FT, 4}
        @test typeof(lb_2) == SMatrix{2, 2, FT, 4}
        @test typeof(ts_2) == SMatrix{3, 2, FT, 6}
        @test eltype(ub_2) == FT
        @test eltype(lb_2) == FT
        @test eltype(ts_2) == FT

        testfunc2(x, y) = x[1] + y[1]
        ub = CNM.UpperOnly{typeof(testfunc2)}(testfunc2)
        lb = CNM.LowerOnly{typeof(testfunc2)}(testfunc2)
        ts = CNM.TwoSided{typeof(testfunc2), typeof(testfunc2)}(
            testfunc2,
            testfunc2,
        )
        inp_3 = SVector{1, FT}(3)
        ans_3 = SMatrix{2, 1, FT}([6; 3])
        ans_33 = SMatrix{3, 1, FT}([6; 6; 3])
        ub_3 = CNM.boundary_connection(ub, inp_3, inp_3)
        lb_3 = CNM.boundary_connection(lb, inp_3, inp_3)
        ts_3 = CNM.boundary_connection(ts, inp_3, inp_3)

        @test ub_3 == ans_3
        @test lb_3 == ans_3
        @test ts_3 == ans_33
        @test typeof(ub_3) == SMatrix{2, 1, FT, 2}
        @test typeof(lb_3) == SMatrix{2, 1, FT, 2}
        @test typeof(ts_3) == SMatrix{3, 1, FT, 3}
        @test eltype(ub_3) == FT
        @test eltype(lb_3) == FT
        @test eltype(ts_3) == FT

        test_matrix = Matrix{FT}([1 -1])
        ml = CNM.MulLayer(test_matrix)
        @test eltype(ml.W) == FT
        @test typeof(ml) === CNM.MulLayer{FT, Matrix{FT}}
        @test ml.W == test_matrix
        @test ml(FT.([-1; 1]))[1] == -2
        @test ml(SMatrix{2, 1, FT}([-1; 2]))[1] == -3

        ml_conv = CNM.convert_model(ml, Float64)
        @test eltype(ml_conv.W) == Float64
        @test typeof(ml_conv) == CNM.MulLayer{Float64, Matrix{Float64}}
        @test ml_conv.W == test_matrix
        @test ml_conv(Float64.([-1; 1]))[1] == -2
        @test ml_conv(SMatrix{2, 1, Float64}([-1; 2]))[1] == -3

        gm_1 = CNM.get_scaling_matrix(ub, FT(-1))
        gm_2 = CNM.get_scaling_matrix(lb, Float64(-1))
        gm_3 = CNM.get_scaling_matrix(ts, FT(-1))
        @test gm_1 == [1 0; 0 -1]
        @test eltype(gm_1) == FT
        @test gm_2 == [1 0; 0 -1]
        @test eltype(gm_2) == Float64
        @test gm_3 == [1 0 0; 0 1 0; 0 0 -1]
        @test eltype(gm_3) == FT

        def_l = Matrix{FT}([1 0; -1 0; -1 1])
        def_u = Matrix{FT}([1 0; -1 0; 1 -1])
        def_t = Matrix{FT}([1 0 0; -1 0 0; 0 1 0; 0 -1 0; 1 0 -1])
        inp_u = Matrix{FT}([1 2 3; 0 4 4])
        inp_l = Matrix{FT}([0 4 4; 1 2 3])
        inp_t = Matrix{FT}([1 1 1; 0 0 0; 2 0.5 -1])
        gd_u = CNM.default_fixed_layers(ub, FT)
        gd_l = CNM.default_fixed_layers(lb, Float64)
        gd_t = CNM.default_fixed_layers(ts, FT)
        @test gd_u[1].weight == def_u
        @test gd_l[1].weight == def_l
        @test gd_t[1].weight == def_t
        @test eltype(gd_u[1].weight) == FT
        @test eltype(gd_l[1].weight) == Float64
        @test eltype(gd_t[1].weight) == FT
        @test gd_u(inp_u) == [0 2 3] #applies upper bound
        @test gd_l(inp_l) == [1 4 4] #applies lower bound
        @test gd_t(inp_t) == [1 0.5 0] #applies both
        test_dec_doc = "this is a test declarative doc\n"
        test_m_doc = "this is a test method doc\n"

        prev_length = length(CNM._BOUND_TYPES_)

        """
        this is a test declarative doc
        """
        @bound_type struct up_functor
            a::Int
            b::DataType
            c::FT
        end

        @bound function (b::up_functor)(
            pred::AbstractArray{T},
            input::AbstractArray{T},
        )::AbstractArray{T} where {T <: AbstractFloat}
            return pred .+ input
        end

        upper_funct = up_functor(1, Tuple{FT(2), FT(3)}, FT(4))

        test_code = "@bound_type struct up_functor\n    a::Int\n    b::DataType\n    c::FT\nend"
        funct_code = "@bound function ((b::up_functor)(pred::AbstractArray{T}, input::AbstractArray{T})::AbstractArray{T}) where T <: AbstractFloat\n    return pred .+ input\nend"
        @test length(CNM._BOUND_TYPES_) == prev_length + 1
        @test CNM.is_bound_type(:up_functor)
        @test CNM.is_bound_type(up_functor)
        @test !CNM.is_bound_type(:low_functor) #not made
        new_info = CNM._BOUND_TYPES_[:up_functor]
        @test !haskey(new_info, :methods)
        CNM.populate_functor_methods!(up_functor)
        @test haskey(new_info, :methods)
        @test length(new_info[:methods]) == 3 #2 constructors (default plus one above), 1 instancemethod
        @test CNM.get_bound_type_info(up_functor) == new_info
        @test new_info[:doc] == test_dec_doc
        @test new_info[:supertype] == Any
        @test new_info[:typevars] == (Int, DataType, FT)
        @test new_info[:code] == test_code
        @test getindex.(new_info[:params], 1) == [:a, :b, :c]
        @test CNM.bound_symbol(up_functor) == :up_functor
        @test CNM.bound_symbol(upper_funct) == :up_functor
        @test CNM.is_valid_bound(:up_functor)
        @test CNM.is_valid_bound(up_functor)
        @test CNM.is_valid_bound(upper_funct)

        m_up = CNM.instancemethods(up_functor)
        @test length(m_up) == 1
        m_up = m_up[1]
        @test m_up in keys(new_info[:methods])
        @test CNM.bound_symbol(upper_funct) == :up_functor
        @test CNM.bound_symbol(:upper_funct) == :upper_funct
        @test !CNM.is_valid_bound(:upper_funct) #can't pass an instance of a callable type as a symbol.
        @test CNM.is_valid_bound(m_up)
        @test CNM.get_bound_info(up_functor) ==
              CNM._BOUND_INFO_[:FUNCTOR][:up_functor]
        up_info = CNM.get_bound_info(upper_funct)
        m_up_info = first(values(CNM._BOUND_INFO_[:FUNCTOR][:up_functor]))
        @test CNM.get_bound_info(m_up) == m_up_info
        @test m_up_info[:method] == m_up
        @test m_up_info[:class] == :FUNCTOR
        @test m_up_info[:type] == (:generic, :dynamic)
        @test isnothing(m_up_info[:docs])
        @test m_up_info[:code] == funct_code
        @test m_up_info[:argtypes] == (
            input = AbstractArray{<:AbstractFloat},
            pred = AbstractArray{<:AbstractFloat},
        )
        @test m_up_info[:ret_type] == AbstractArray{<:AbstractFloat}

        prev_func_length = length(CNM._BOUND_INFO_[:FUNCTION])

        """
        this is a test method doc
        """
        @bound function lower(
            pred::AbstractArray{<:AbstractFloat},
            input::AbstractArray{<:AbstractFloat},
        )::AbstractArray{<:AbstractFloat}
            return pred .- input
        end
        func_code = "@bound function lower(pred::AbstractArray{<:AbstractFloat}, input::AbstractArray{<:AbstractFloat})::AbstractArray{<:AbstractFloat}\n    return pred .- input\nend"
        @test length(CNM._BOUND_INFO_[:FUNCTION]) == prev_func_length + 1
        @test CNM.bound_symbol(lower) == :lower
        m_low = methods(lower)
        @test length(m_low) == 1
        m_low = m_low[1]
        @test CNM.bound_symbol(m_low) == :lower
        @test CNM.bound_symbol(:lower) == :lower
        @test CNM.is_valid_bound(lower)
        @test CNM.is_valid_bound(m_low)
        @test CNM.is_valid_bound(:lower)
        @test !CNM.is_bound_type(:lower) #it is not a functor
        @test CNM.get_bound_info(lower) == CNM._BOUND_INFO_[:FUNCTION][:lower]
        low_info = CNM.get_bound_info(:lower)
        m_low_info = first(values(CNM._BOUND_INFO_[:FUNCTION][:lower]))
        @test CNM.get_bound_info(m_low) == m_low_info
        @test m_low_info[:method] == m_low
        @test m_low_info[:class] == :FUNCTION
        @test m_low_info[:type] == (:generic, :dynamic)
        @test m_low_info[:docs] == test_m_doc
        @test m_low_info[:code] == func_code
        @test m_low_info[:argtypes] == (
            input = AbstractArray{<:AbstractFloat},
            pred = AbstractArray{<:AbstractFloat},
        )
        @test m_low_info[:ret_type] == AbstractArray{<:AbstractFloat}

        inps, classes = CNM.get_bound_evaluation_modes(low_info)
        @test inps == [:generic]
        @test classes == [:dynamic]
        @test isnothing(CNM.check_evaluation_mode(inps, m_low)) #generic bound should yield nothing
        inps, classes = CNM.get_bound_evaluation_modes(up_info)
        @test inps == [:generic]
        @test classes == [:dynamic]
        @test isnothing(CNM.check_evaluation_mode(inps, m_up)) #generic bound should yield nothing

        up_const = CNM.buildUpper(upper_funct)
        @test typeof(up_const) <: CNM.UpperOnly
        @test up_const.bound == upper_funct
        low_const = CNM.buildLower(lower)
        @test typeof(low_const) <: CNM.LowerOnly
        @test low_const.bound == lower
        both_const = CNM.buildTwoSided(upper_funct, lower)
        @test typeof(both_const) <: CNM.TwoSided
        @test both_const.upper_bound == upper_funct
        @test both_const.lower_bound == lower

        up_const2 = CNM.buildConstraints(upper_funct, nothing)
        @test up_const2 == up_const
        low_const2 = CNM.buildConstraints(nothing, lower)
        @test low_const2 == low_const
        both_const2 = CNM.buildConstraints(upper_funct, lower)
        @test both_const2 == both_const

        @test CNM.has_evaluation_mode(up_const2, :dynamic)
        @test CNM.has_evaluation_mode(low_const2, :dynamic)
        @test CNM.has_evaluation_mode(both_const2, :dynamic)
        @test !CNM.has_evaluation_mode(up_const2, :static)
        @test !CNM.has_evaluation_mode(low_const2, :static)
        @test !CNM.has_evaluation_mode(both_const2, :static)

        struct TestType
            f1::Int
            f2::FT
        end
        test_data = Dict(:name => :TestType, :code => "no code")
        CNM.construct_type_info!(TestType, test_data)
        test_info = CNM._BOUND_TYPES_[:TestType]
        @test test_info[:code] == "no code"
        Base.delete!(CNM._BOUND_TYPES_, :TestType)

        @test CNM.get_methods(up_functor) == CNM.instancemethods(up_functor)
        @test CNM.get_methods(lower) == methods(lower)

        test_string_1 = "yip#= REPL[210]:1 =#pee"
        test_string_2 = """yip#= c:\\Users\\bob\\Documents\\Caltech\\Research\\ClimaLandTesting.jl:10 =#pee"""
        @test CNM._strip_location_comments(test_string_1) == "yippee"
        @test CNM._strip_location_comments(test_string_2) == "yippee"

        f1 = :(function f1(x::DataType, y::Array{N})::G where {G, N, T <: Tuple}
            return x * y
        end)
        f2 = :(g(y, z::Float32) = 3)
        f3 = :(
            function (::MyType{Q, P})(
                a::Tuple{Tuple{G, N}, Array{K}},
                b::Array,
            )::A where {A <: AbstractFloat, G, N}
                return a
            end
        )
        n1, t1, out_1 = CNM._get_func_info(f1)
        n2, t2, out_2 = CNM._get_func_info(f2)
        n3, t3, out_3 = CNM._get_func_info(f3)
        @test [n1, n2, n3] == [:f1, :g, :MyType]
        @test [t1, t2, t3] == [:FUNCTION, :FUNCTION, :FUNCTOR]
        @test eval(out_1) == Any
        @test isnothing(eval(out_2))
        @test eval(out_3) == AbstractFloat

        """
        I am making a variable called x
        """
        test_declare_var = 3

        #Would probably appreciate some help in coming up with more pathological examples,
        #to make sure this function is accomplishing what we want it to:
        @test CNM._get_argtypes(m_low.sig) == (
            typeof(lower),
            AbstractArray{<:AbstractFloat},
            AbstractArray{<:AbstractFloat},
        )
        #parametric types without explicitly defined T in codespace means we need to compare these as strings,
        #but this tests the "transform" capability as well of the function:
        @test CNM._get_argtypes(m_up.sig, string) == (
            "up_functor",
            "AbstractArray{T<:AbstractFloat}",
            "AbstractArray{T<:AbstractFloat}",
        )
        @test CNM._get_argtypes(Union{}) == ()
        @test CNM._get_argtypes(up_functor) == fieldtypes(up_functor)
        @test CNM._get_declarative_docs(up_functor) == test_dec_doc
        @test CNM._get_declarative_docs(:test_declare_var) ==
              "I am making a variable called x\n"

        @test CNM._get_doc_from_method(m_low) == test_m_doc
        @test isnothing(CNM._get_doc_from_method(m_up))

        mt = SMatrix{1, N, <:AbstractFloat, N} where {N <: Int}

        st_expr = :(function static(pred::mt, input::mt)::mt
            return pred .* inp
        end)
        st_name, st_type, st_ret_expr = CNM._get_func_info(st_expr)
        eval(st_expr) #make the function
        static_m = methods(static)[1]

        st_data = Dict(
            :code => "no code again",
            :name => st_name,
            :class => st_type,
            :ret_type => eval(st_ret_expr),
        )

        CNM.construct_method_info!(static_m.name, static_m, st_data)
        @test CNM.is_valid_bound(static_m)
        grab_static = CNM.get_bound_info(static_m)
        @test grab_static[:ret_type] == mt
        @test grab_static[:code] == "no code again"
        @test grab_static[:name] == :static
        @test grab_static[:class] == :FUNCTION
        @test grab_static[:type] == (:batched, :static) #validates eval_modes funcs again too
        @test isnothing(grab_static[:docs])
        @test grab_static[:argtypes].pred == mt
        @test grab_static[:argtypes].input == mt
        Base.delete!(CNM._BOUND_INFO_[:FUNCTION], :static)

        @test CNM.get_val(Val{Float32}()) == FT
        @test CNM.get_val(Val{true}())
        @test CNM.get_val(Val{3}()) == 3

        #TODO:
        #ConstrainedNeuralModel type
        #ConstrainedNeuralModel constructor
        #ConstrainedNeuralModel Functor methods, both types
        #Flux.trainable() on different CNMs
        #set_predictive_model_out_scale!
        #scale_model!
        #save_model (how to test this one?)
        #load_model in various forms (requires the NeuralSnow artifacts to be saved)
        #inspect_model_metadata
        #make_static_model
        #make_mutable_dynamic_model
        #make_mutable_fixed_model
        #trainmodel

        #get_modules
        #_check_children
        #build_parameter_metadata
        #assess_model_transferability
        #build_API
        #get_fixed_layer_info
        #build_bound_docs
        #build_model_metadata

        #clear the test types from the _BOUND_INFO_, _BOUND_TYPES_ when you are done testing
        #Delete/fix relevant NeuralSnow tests
    end

    @testset "Testing NeuralSnow module" begin
        #Model setup:
        FT = Float32
        toml_dict = LP.create_toml_dict(FT)
        start_date = DateTime(2005)
        Δt = FT(180.0)
        domain = Point(; z_sfc = FT(0))
        "Radiation"
        SW_d = TimeVaryingInput((t) -> eltype(t)(20.0))
        LW_d = TimeVaryingInput((t) -> eltype(t)(20.0))
        rad = ClimaLand.PrescribedRadiativeFluxes(FT, SW_d, LW_d, start_date)
        "Atmos"
        precip = TimeVaryingInput((t) -> eltype(t)(0.01 / Δt))
        T_atmos = TimeVaryingInput((t) -> eltype(t)(278.0))
        u_atmos = TimeVaryingInput((t) -> eltype(t)(10.0))
        q_atmos = TimeVaryingInput((t) -> eltype(t)(0.03))
        h_atmos = FT(3)
        P_atmos = TimeVaryingInput((t) -> eltype(t)(101325))
        atmos = ClimaLand.PrescribedAtmosphere(
            precip,
            precip,
            T_atmos,
            u_atmos,
            q_atmos,
            P_atmos,
            start_date,
            h_atmos,
            toml_dict,
        )

        #Test extension utilities
        z_model = NeuralSnow.get_znetwork(FT)
        @test typeof(z_model) <: Flux.Chain
        z32 = NeuralSnow.converted_model_type(z_model, FT)
        @test eltype(z32[1].layers[1].weight) == FT
        z64 = NeuralSnow.converted_model_type(z_model, Float64)
        @test eltype(z64[1].layers[1].weight) == Float64
        dens_model1 = NeuralSnow.NeuralDepthModel(FT)
        @test prod(
            Flux.trainables(dens_model1.z_model) .== Flux.trainables(z32),
        )
        @test eltype(dens_model1.α) == FT
        test_alph = 3 / 86400
        dens_model2 = NeuralSnow.NeuralDepthModel(
            FT,
            α = test_alph,
            Δt = Δt;
            model = NeuralSnow.get_znetwork(FT),
        )
        @test dens_model2.α == FT(test_alph)
        @test dens_model2.z_model[:final_scale].weight[2, 2] == FT(1 / Δt)

        parameters = SnowParameters(toml_dict, Δt; density = dens_model2)
        model = ClimaLand.Snow.SnowModel(
            parameters = parameters,
            domain = domain,
            boundary_conditions = ClimaLand.Snow.AtmosDrivenSnowBC(atmos, rad),
        )

        drivers = ClimaLand.get_drivers(model)
        Y, p, coords = ClimaLand.initialize(model)
        @test (Y.snow |> propertynames) ==
              (:S, :S_l, :U, :Z, :P_avg, :T_avg, :R_avg, :Qrel_avg, :u_avg)

        Y.snow.S .= FT(0.1)
        Y.snow.U .=
            ClimaLand.Snow.energy_from_T_and_swe.(
                Y.snow.S,
                FT(273.0),
                Ref(model.parameters),
            )
        Y.snow.Z .= FT(0.2)
        set_initial_cache! = ClimaLand.make_set_initial_cache(model)
        t0 = FT(0.0)
        set_initial_cache!(p, Y, t0)
        oldρ = p.snow.ρ_snow
        NeuralSnow.update_density_and_depth!(
            p.snow.ρ_snow,
            p.snow.z_snow,
            model.parameters.density,
            Y,
            p,
            model.parameters,
        )
        @test p.snow.z_snow == Y.snow.Z
        @test p.snow.ρ_snow == oldρ
        output1 = NeuralSnow.eval_nn(dens_model2, FT.([0, 0, 0, 0, 0, 0, 0])...)

        @test eltype(output1) == FT
        @test output1 == 0.0f0

        zerofield = similar(Y.snow.Z)
        zerofield .= FT(0)
        dY = similar(Y)
        NeuralSnow.update_dzdt!(dY.snow.Z, dens_model2, Y)
        @test dY.snow.Z == zerofield

        Z = FT(0.5)
        S = FT(0.1)
        dzdt = FT(1 / Δt)
        dsdt = FT(1 / Δt)
        @test NeuralSnow.clip_dZdt(S, Z, dsdt, dzdt, Δt) == dzdt

        @test NeuralSnow.clip_dZdt(Z, S, dsdt, dzdt, Δt) ≈ FT(1.4 / Δt)

        @test NeuralSnow.clip_dZdt(S, Z, FT(-S / Δt), dzdt, Δt) ≈ FT(-Z / Δt)


        dswe_by_precip = 0.1
        Y.snow.P_avg .= FT(dswe_by_precip / Δt)
        exp_tendency! = ClimaLand.make_compute_exp_tendency(model)
        exp_tendency!(dY, Y, p, FT(0.0))
        @test Array(parent(dY.snow.Z))[1] * Δt > dswe_by_precip
        new_dYP = FT(test_alph) .* (p.drivers.P_snow .- Y.snow.P_avg)
        @test dY.snow.P_avg == new_dYP

    end
end
