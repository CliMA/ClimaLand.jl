# Solver Functions

```@meta
CurrentModule = ClimaLand
```

## Prognostic Variables

```@docs
ClimaLand.prognostic_vars
ClimaLand.prognostic_types
ClimaLand.prognostic_domain_names
ClimaLand.initialize_prognostic
```

## Sources

```@docs
ClimaLand.AbstractSource
ClimaLand.source!
```

## Boundary Conditions

```@docs
ClimaLand.AbstractBC
ClimaLand.AbstractBoundary
ClimaLand.TopBoundary
ClimaLand.BottomBoundary
ClimaLand.boundary_vars
ClimaLand.boundary_var_domain_names
ClimaLand.boundary_var_types
```

## Implicit Tendencies

```@docs
ClimaLand.make_imp_tendency
ClimaLand.make_compute_imp_tendency
ClimaLand.make_compute_jacobian
ClimaLand.make_jacobian
ClimaLand.FieldMatrixWithSolver
```

## Explicit Tendencies

```@docs
ClimaLand.make_exp_tendency
ClimaLand.make_compute_exp_tendency
```
