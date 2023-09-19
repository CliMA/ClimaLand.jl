#=
# Domain Tutorial

## Goals of the tutorial
The goal of this is to outline what is currently implemented in ClimaLSM
and to serve as a software design document
for future development involving the underlying domains.

## Background

In both the atmosphere and the ocean, all variables are defined at all locations
in the region of interest, or domain.  For example, the air density, temperature, pressure,
and wind speed are defined everywhere in the domain. After choosing a resolution
and discretizing space, the numerical problem is to advance a
system of differential equations, where at each coordinate point a value of
`ρ`, `T`, `P`, and `u⃗` are solved for at each step. The choice of domain is a question "only"
of geometry: you may be interested in a large eddy simulation (using a box domain), or
in a global model (where you would need a spherical shell domain
representing the atmosphere or ocean from some depth to z_sfc = 0).

For land surface
models, each variable is not defined everywhere in space. For example,
the soil water content `θ` is only defined below ground. Snow water equivalent (`S`) is only
defined on the surface itself. Canopy variables are only defined above ground.
Once we have discretized the land surface region into
a set of points, the numerical problem is to advance a system of ODEs, where
at each coordinate point a different subset of (`θ`, `S`, ...) are solved for.

In other words, different variables in land surface models exist
in different, overlapping, domains. We need to decide on the geometry of interest (e.g. single column
vs a global simulation), but we also need to specify where each variable of the model is defined.

ClimaLSM Domains were designed with this in mind. The domains are defined
so that
1. the user can easily switch geometries, e.g. single column to global model,
2. individual component models can be run by themselves, using a single domain,
3. the same domains can be used to set up multi-component models (LSMs),
4. different variables can exist on different parts of the domain.


## What is a ClimaLSM Domain?
A domain represents a region of space. In ClimaLSM, domains are simply
structs containing parameters that define these regions - for example an
x-range and y-range that define a plane. In addition, ClimaLSM domains store
the ClimaCore function spaces for the
physical domain as a NamedTuple. When solving partial differential equations, the spatial
discretization is tied to a set of basis functions you wish to use to represent the prognostic
variable as a function of space. The nodal points - the locations in space where the variable
is solved for - are arranged in `space` in a manner which depends on these basis functions.
Note that these spaces are only mathematically needed when your variables satisfy PDEs[^1],
but that they still exist when your variables do not, because we are using the same 
underlying infrastructure in both cases.


## Domain types
All ClimaLSM domains are subtypes of abstract type
[`ClimaLSM.Domains.AbstractDomain`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.Domains.AbstractDomain).
A variety of concrete domain types are supported:

- 0D: [`Domains.Point`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.Domains.Point)
- 1D: [`Domains.Column`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.Domains.Column)
- 2D: [`Domains.Plane`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.Domains.Plane), [`Domains.SphericalSurface`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.Domains.SphericalSurface)
- 3D: [`Domains.HybridBox`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.Domains.HybridBox), [`Domains.SphericalShell`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.Domains.SphericalShell).

As discussed above, our modeling requires that variables of a model can be defined on different subsets of the domain. Because of that, we define the concept of a surface domain, and a subsurface domain. Not all domains have a surface and subsurface; some only have surface domains, as shown in the Table below.


|  Domain | Surface Domain | Subsurface Domain |
| :---         |     :---:      |          ---: |
| Column | Point    | Column    |
| HybridBox    | Plane      | HybridBox      |
| SphericalShell    | SphericalSurface      | SphericalShell      |


There is a single key method which take a ClimaLSM domain as an argument.

- [` coordinates(domain)`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.Domains.coordinates): under the hood, this function  uses
the NamedTuple of function spaces (domain.space) to create the coordinate field for the surface and subsurface domains (as applicable), stored in a NamedTuple.
Depending on the domain, the returned coordinate field will have elements of different names and types. For example,
the SphericalShell domain has subsurface coordinates of latitude, longitude, and depth, while the surface coordinates
are latitude and longitude. A Plane domain has coordinates
of x and y (surface only), and a Point domain only has a coordinate z_sfc (surface only). Column domains have a surface coordinate of z_sfc,
and subsurface coordinates of z.

It is important to note that the horizontal domain used for the surface and subsurface
domains are identical in all simulations. This ensures that we can use the same indexing
of surface and subsurface domains and variables. Otherwise we would need
to develop additional infrastructure in order to, for example, select the correct subsurface
column corresponding to a particular surface location.


## How variable initialization depends on domains
Single component models (soil, snow, vegetation, canopy...) must have an associated domain in order
to solve the their equations.  Which domain is appropriate depends on the model equations and
on the configuration of interest (single column or global, etc.). For example, the soil model is a
vertically resolved model, so only domains with a vertical extent (Column, HybridBox, or SphericalShell)
make sense to use. A single layer snow model does not require vertical resolution - and so the domains
that make sense to use are a Point, Plane, or SphericalSurface.

When a developer first defines a model, they need to specify the symbols used for the prognostic variables,
via [`prognostic_vars`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.prognostic_vars),
and
the types of those variables,
 via [`prognostic_types`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.prognostic_types).

They additionally need to define which subset of the domain the variables are defined on, using 
[`prognostic_domain_names`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.prognostic_domain_names).

The [`initialize`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.initialize)
function (which calls both
[`initialize_prognostic`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.initialize_prognostic)
 and [`initialize_auxiliary`](https://clima.github.io/ClimaLSM.jl/dev/APIs/shared_utilities/#ClimaLSM.initialize_auxiliary))
creates the prognostic state vector `Y` (a ClimaCore.Fields.FieldVector). Each field (ClimaCore.Fields.Field) stored within
the field vector corresponds to a prognostic variable (identified with the symbol specified). If the prognostic type for that variable
is a float, the field will be a field of float values (a scalar field)[^4].

How do domains tie into this? The field of a prognostic variable corresponds in a 1-1 fashion with the coordinate field of the subset of the domain associated with that variable via `prognostic_domain_name`.  For example, the bucket model has a vertically resolved temperature `T`, but the bucket water content 
`W` is not vertically resolved. If your domain is a Column, the subsurface coordinates may be [-4.5,-3.5,-2.5,-1.5, -0.5], and 
the surface coordinate would be [-0.0]. Your prognostic variable field for `T` will be [T[-4.5], T[-3.5]; T[-2.5], T[-1.5], T[-0.5]], and for `W`  it will be [W[0.0],]. Your variable always has the same spatial resolution as the associated subset of the domain. 

This functionality is not required for every standalone component model. For example, a single layer snow model
will only have variables on the surface of the domain (which in this case, would be the entire Point, Plane, or 
SphericalShell domain). The user still must define the prognostic_domain_names method. This functionality is required
for most multi-component models.


## Future work
Almost all interactions between variables in land surface models are within column - that is, there is only
vertical transport and exchanges. The exception to this is the horizontal flow of water on the surface
and within the soil. The tendency (produced by `make_exp_tendency` and `make_imp_tendency`) functions
(the ODE functions) can be split into "vertical" and "horizontal" pieces.

We envision each step of the land surface model simulation to be solved  in two steps: (1) the vertical tendency
evaluations are carried out (and can be parallelized), and (2) the horizontal tendency functions are then evaluated
(possibly less frequently?) and require communcation between columns.
In this case, tendency functions will need to be aware of the domain.
In general, tendencies reflecting horizontal flow will be treated explicitly and include in the explicit tendency function.
Tendencies reflecting vertical flow may be treated explicitly or implicitly depending on the use case. To solve the problem,
we then use IMEX (mixed explicit/implicit) methods.


[^1]: finite differencing is used in the vertical, and spectral elements are used in the horizontal.

[^2]: a suprasurface region may also be necessary - for example if the canopy airspace model involves PDEs.

[^3]: We also will support having an array-like type of variable.
=#
