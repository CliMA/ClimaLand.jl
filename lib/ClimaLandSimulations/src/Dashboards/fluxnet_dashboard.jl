export fluxnet_dashboard, fluxnet_app

"""
    fluxnet_dashboard(menu, menu2)

Return selected figure from selected site. 
"""
function fluxnet_dashboard(menu, menu2)
    display_fig = map(menu.value) do value
        @info value
        sv = run_fluxnet(value)[1]
        simulation_inputs = make_inputs_df(value)[1]
        simulation_outputs = make_output_df(value, sv, simulation_inputs)
        figs = make_plots(
            simulation_inputs,
            simulation_outputs;
            save_fig = false,
            dashboard = true,
        )
        return Dict(
            "timeseries" => figs.timeseries,
            "water" => figs.water,
            "fingerprint" => figs.fingerprint,
            "diurnals" => figs.diurnals,
        )
    end

    fig = map(menu2.value) do value
        @info value
        return get(display_fig[], value, nothing)
    end
    return fig
end

"""
    fluxnet_app()

Make a web dashboard to interact with fluxnet_dashboard(menu, menu2)
"""
function fluxnet_app()
    ClimaLand_dashboard = App() do
        sites_ID = ["US-MOz", "US-Ha1", "US-NR1", "US-Var"]
        menu = Dropdown(sites_ID)
        displayed_figure =
            ["timeseries", "water", "fingerprint", "diurnals"]
        menu2 = Dropdown(displayed_figure)
        fig = fluxnet_dashboard(menu, menu2)
        return DOM.div(menu, menu2, fig)
    end
    return ClimaLand_dashboard
end
