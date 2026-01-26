# Public APIs

## `ClimaDiagnostics`

```@docs
ClimaDiagnostics.DiagnosticsHandler
ClimaDiagnostics.orchestrate_diagnostics
ClimaDiagnostics.DiagnosticsCallback
ClimaDiagnostics.IntegratorWithDiagnostics
```

## `Schedules`

```@docs
ClimaDiagnostics.Schedules.AbstractSchedule
ClimaDiagnostics.Schedules.short_name
ClimaDiagnostics.Schedules.long_name
ClimaDiagnostics.Schedules.DivisorSchedule
ClimaDiagnostics.Schedules.EveryStepSchedule
ClimaDiagnostics.Schedules.EveryDtSchedule
ClimaDiagnostics.Schedules.EveryCalendarDtSchedule
ClimaDiagnostics.Schedules.EveryCalendarDtSchedule(dt, t::ClimaUtilities.TimeManager.ITime)
```

## `DiagnosticVariables`

```@docs
ClimaDiagnostics.DiagnosticVariables.DiagnosticVariable
ClimaDiagnostics.DiagnosticVariables.short_name
ClimaDiagnostics.DiagnosticVariables.long_name(dv::ClimaDiagnostics.DiagnosticVariables.DiagnosticVariable)
ClimaDiagnostics.DiagnosticVariables.descriptive_short_name
ClimaDiagnostics.DiagnosticVariables.descriptive_long_name
```

## `ScheduledDiagnostics`

```@docs
ClimaDiagnostics.ScheduledDiagnostics.ScheduledDiagnostic
ClimaDiagnostics.ScheduledDiagnostics.output_short_name
ClimaDiagnostics.ScheduledDiagnostics.output_long_name
```


## `Writers`

```@docs
ClimaDiagnostics.AbstractWriter
ClimaDiagnostics.Writers.DictWriter
ClimaDiagnostics.Writers.NetCDFWriter
ClimaDiagnostics.Writers.HDF5Writer
ClimaDiagnostics.Writers.interpolate_field!
ClimaDiagnostics.Writers.write_field!
ClimaDiagnostics.Writers.default_num_points
Base.close
```
