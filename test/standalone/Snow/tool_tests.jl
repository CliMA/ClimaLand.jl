using ClimaLSM.Snow.DataTools
using ClimaLSM.Snow.ModelTools
using Test
using BSON, Dates, HTTP
using DataFrames, CSV, StatsBase, Flux, LinearAlgebra

@testset "Testing Data Utilities" begin
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

    test_data2 = sitedata_daily(
        station_id,
        station_state,
        start = start_date,
        finish = end_date,
    )
    @test typeof(test_data2) === DataFrame
    @test size(test_data2) == (32, 8)
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
    @test DataFrames.names(test_data2) == colnames
    @test typeof(test_data2[1, 1]) == Date
    @test sum(test_data2[!, :z]) == 1168

    download_link = "https://caltech.box.com/shared/static/4tih9hiydrc7bcvrpkr2x727b5l54oiq.csv"
    download_data = CSV.read(HTTP.get(download_link).body, DataFrame)
    @test isequal(download_data, test_data2)

    test_data3 = sitedata_hourly(
        station_id,
        station_state,
        start = start_date,
        finish = end_date,
    )
    @test typeof(test_data3) === DataFrame
    @test size(test_data3) == (768, 8)
    @test DataFrames.names(test_data3) == colnames
    @test typeof(test_data3[1, 1]) == DateTime
    @test sum(skipmissing(test_data3[!, :z])) == 25063
    @test sum(skipmissing(test_data3[!, :rel_hum_avg])) == 62232

    test_bounds = Dict{Symbol, Tuple{Real, Real}}(
        :z => (0, 37),
        :rel_hum_avg => (84, 100),
    )
    test_data4 = apply_bounds(test_data3, test_bounds)
    @test sum(skipmissing(test_data4[!, :z])) == 19007
    @test sum(skipmissing(test_data4[!, :rel_hum_avg])) == 45973

    test_data5 = hourly2daily(test_data3)
    @test sum(test_data5[!, :z]) == 1168
    @test sum(test_data5[!, :rel_hum_avg]) == 2593
    @test sum(describe(test_data5, :nmissing)[!, 2]) == 0

    test_data6 = rectify_daily_hourly(test_data2, test_data5)
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
    test_data7 = scale_cols(test_data6, scales)
    @test sum(test_data7[!, :rel_hum_avg]) ≈ 25.93 atol = 1e-3
    test_data8 = makediffs(test_data7, Day(1))
    @test size(test_data8) == (31, 11)
    @test mean(test_data8[!, :dzdt]) ≈ 0 atol = 1e-7
    @test mean(test_data8[!, :dSWEdt]) ≈ 0 atol = 1e-7
    @test minimum(test_data8[!, :dprecipdt]) == 0

    test_data9 = rolldata(test_data8, Day(1), 7)
    @test size(test_data9) == (25, 11)
    @test mean(test_data9[!, :dzdt]) ≈ 0 atol = 1e-8
    @test mean(test_data9[!, :dSWEdt]) ≈ 0 atol = 1e-8
    @test minimum(test_data9[!, :dprecipdt]) == 0
    @test minimum(test_data9[!, :z]) ≈ 0.8636 atol = 1e-4

    test_data8[!, :id] .= station_id
    test_data10 =
        prep_data(test_data8, extract_vars = [:dprecipdt_rain, :dprecipdt_snow])
    x_train, y_train =
        make_data(test_data10, [:dprecipdt_rain], :dprecipdt_snow, 1.0)
    @test test_data10[!, 1] .+ test_data10[!, 2] == test_data8[!, :dprecipdt]
    @test size(x_train) == (1, 31)
    @test size(y_train) == (1, 31)
end

@testset "Testing Model Utilities" begin
    nfeatures = 7
    n = 5
    pred_vars = [
        :z,
        :SWE,
        :rel_hum_avg,
        :sol_rad_avg,
        :wind_speed_avg,
        :air_temp_avg,
        :dprecipdt_snow,
    ]
    z_idx = 1
    p_idx = 7
    model = make_model(nfeatures, n, z_idx, p_idx)
    ps = get_model_ps(model)
    for item in ps
        item[:] .= Float32(1.0)
    end
    test_input = Matrix{Float32}(ones(nfeatures, 8))
    @test model(test_input)[1] == 1968
    @test size(model(test_input)) == (1, 8)
    @test sum(length, ps) == nfeatures * (n * (2 * nfeatures + 1) + 2) + 1
    setoutscale!(model, 0.5)
    @test model[:final_scale].weight[3, 3] == 0.5
    @test model(test_input)[1] == 984
    setoutscale!(model, 2.0)
    @test model(test_input)[1] == 1968 #this tests clipping
    @test model(-test_input)[1] == 0.0
    settimescale!(model, 1968)
    @test model[:final_scale].weight[2, 2] == Float32(1 / 1968)
    @test ModelTools.evaluate(model, Vector{Float32}(ones(nfeatures)))[1] ==
          1968

    test_constants = [1.0, 2.0, 3.0, 4.0, 5.0]
    x = zeros(6, 5)
    x[1:5, 1:5] = diagm([1, 2, 3, 4, 5])
    x[6, :] = [1, 1, 1, 1, 1]
    y = [1.0, 4.0, 9.0, 16.0, 25.0, 15.0]
    x_dataframe = DataFrame(x, :auto)
    x_dataframe[!, :y] = y
    answer = [1.0, 2.0, 3.0, 4.0, 5.0, 0.0]
    model2 = LinearModel(x_dataframe, [:x1, :x2, :x3, :x4, :x5], :y)
    model3 = LinearModel(x, y)
    @test model2 ≈ answer
    @test model3 ≈ answer
    @test ModelTools.evaluate(model2, Vector{Float32}(ones(5)))[1] == 15
    @test typeof(ModelTools.evaluate(model2, Vector{Float32}(ones(5)))[1]) ==
          Float32

    temp(x) = x
    test_loss_check(x, y) = ModelTools.custom_loss(x, y, temp, 2, 1)
    input_x = Vector{Float32}([1, 2, 3, 4, 5])
    input_y = Vector{Float32}([2, 3, 5, 2, 3])
    @test test_loss_check(input_x, input_y) == 11.8f0
    @test typeof(test_loss_check(input_x, input_y)) == Float32

    data_download_link = "https://caltech.box.com/shared/static/n59m3iqcgr60gllp65rsrd3k0mtnsfmg.csv"
    model_download_link = "https://caltech.box.com/shared/static/bbu12b518i49aj3pl6b8t9t05twl5teq.bson"
    data = CSV.read(HTTP.get(data_download_link).body, DataFrame)
    data = data[data[!, :id] .== 1286, :]
    data = prep_data(data)
    nmodel = make_model(nfeatures, n, z_idx, p_idx)
    model_state =
        BSON.load(IOBuffer(HTTP.get(model_download_link).body))[:model_state]
    Flux.loadmodel!(nmodel, model_state)
    settimescale!(model, 86400.0)
    pred_series, _, _ = make_timeseries(nmodel, data, Day(1))
    true_series = data[!, :z]
    test_loss(x, y) = ModelTools.custom_loss(x, y, nmodel, 2, 1)
    series_err =
        sqrt(sum((pred_series .- true_series) .^ 2) ./ length(pred_series))
    direct_err =
        test_loss(Matrix{Float32}(select(data, pred_vars))', data[!, :dzdt]')
    @test series_err ≈ 0.1 atol = 0.05
    @test direct_err ≈ 0 atol = 1e-12

    out_scale = maximum(abs.(data[!, :dzdt]))
    x_train, y_train = make_data(data, pred_vars, :dzdt, out_scale)
    ps = get_model_ps(nmodel)
    settimescale!(nmodel, 86400 * out_scale)
    setoutscale!(nmodel, 1.0)
    callback_check = [0.0]
    function call_check(val = callback_check)
        val[1] += 1
    end
    nepochs = 10
    trainmodel!(
        nmodel,
        ps,
        x_train,
        y_train,
        2,
        1,
        nepochs = nepochs,
        cb = call_check,
    )
    @test callback_check[1] == nepochs
    setoutscale!(nmodel, out_scale)
    settimescale!(nmodel, 86400)
    pred_series, _, _ = make_timeseries(nmodel, data, Day(1))
    series_err =
        sqrt(sum((pred_series .- true_series) .^ 2) ./ length(pred_series))
    @test series_err <= 0.2
end
