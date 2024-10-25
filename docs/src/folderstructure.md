# ClimaLand folder structure

[ClimaLand](https://github.com/CliMA/ClimaLand.jl) home directory has 5 main folders:

- docs: contains files to generate [the documentation website](https://clima.github.io/ClimaLand.jl/dev/).
- experiments: contains simple runs of `ClimaLand` models. 
- parameters: contains a file to retrieve constants such as avogadro's number, the speed of light, etc. 
- src: contains the code of `ClimaLand` models. 
- test: contains [unit tests](https://en.wikipedia.org/wiki/Unit_testing), which are meant to ensure small pieces of ClimaLand source code work as intended before merging pull requests.

and 3 GitHub actions folders. GitHub actions are `.yml` files, which are bash scripts that runs on a remote computer on each git push. 
- .buildkite: contains a script building outputs such as figures from experiments and tests folders. These runs are carried out as part of CI and must run without error in order to merge a PR. 
- .dev: contains useful tools for developers, such as a format checker for Julia (which is run as part of CI and must pass before a PR can be merged into main). 
- .github: contains various scripts, for example, this documentation website is built each time a change is push to `ClimaLand`. 

as well as 5 files:
- .gitignore: [commonly used git file](https://git-scm.com/docs/gitignore), contains files, files type, and folders that should be ignored by git. 
- LICENSE: License file of `ClimaLand`, you can [read it](https://github.com/CliMA/ClimaLand.jl/blob/main/LICENSE) to learn about legal practice regarding use of `ClimaLand` open source code. 
- Project.toml: The Julia programming language requires a `Project.toml` file to create an environment, which specify dependencies of a project as well as its version, name, authors and a unique identifier number (uuid). Every Julia registered package has a `Project.toml` file. 
- README.md: This markdown file contains the info that you can read on [ClimaLand GitHub web page](https://github.com/CliMA/ClimaLand.jl)


## /docs folder

Julia packages are recommended to have a `\docs` folder that builds a standardised documentation following the official documentation generator for Julia: [Documenter.jl](https://documenter.juliadocs.org/stable/).

The folder /docs contains:
- a `/src` folder: It is recommended to put your markdown pages inside this folder. Each markdown file (.md extension text file) is a page accessible through the menu of the documentation. For example, `docs/src/Contributing.md` contains the text you can read on the documentation "contribution guide" menu. The path to this .md file and the name of the menu is set in the `docs/make.jl` file. 
- a `make.jl` file: This Julia file contains your documentation website structure. Running this file will build your website pages, you can run it locally, but it is commonly built remotely via .github/workflows/docs.yml to generate the github static page hosted on the [gh-pages branch](https://github.com/CliMA/ClimaLand.jl/tree/gh-pages). 

Note: the documentation can have submenu. For example, `APIs` have submenu `ClimaLand` which has many submenu... This structure is built in our current framework via a file `docs/list_of_apis.jl` in that example, which is then included in `docs/make.jl`.   

## /experiments folder

The `experiments` folder contains scripts to run ClimaLand models. It contains a folder for `integrated` models and a folder for `standalone` models. It is meant to provide users with simple examples of ClimaLand runs. The files contains meteorological inputs (such as precipitation), values for every parameters, and the domains and timestepper are specified. 

For example, `/experiments/LSM/ozark/` contains:
- ozark_domain.jl: Describes the soil domain (depth, number of layer), and the canopy (number and height of stems and leaves).
- ozark_met_drivers_FLUXNET.jl: This files load meteorological input data from the ozark FLUXNET file, and does additional things such as spline interpolation of these drivers.
- ozark_parameters.jl: In this file, parameters values are defined. 
- ozark_simulation.jl: In this file, initial and final time are set, as well as time resolution and time stepper algorithm. 
- ozark.jl: running this script will include all the above scripts, and run ClimaLand for the single-site ozark. It will produce output in a text file as well as some figures comparing data and simulation.

## /src folder

The `/src` folder contains the source code of `ClimaLand` models. It contains 3 folders:
- shared_utilities: This is a core folder that defines functions and data structures used across all modules and models types of `ClimaLand`. For example, `shared_utilities/models.jl` defines and export the function `make_update_aux` which will be used to create a function which updates the auxiliary parameters, stored in the vector `p`, `shared_utilities/boundary_conditions.jl` defines functions for setting boundary condition for PDE domains, etc.
- standalone: This folder contains standalone models, which are submodels that can be run independently of each other. This is an important aspect of `ClimaLand` code design: to maximize modularity, sub-models can be run alone, and many different methods of the same sub-model can be defined via Julia multiple dispatch. The standalone folder is independent from the integrated folder. 
- integrated: This folder contains integrated models. It assembles standalone models together, as one would assemble pieces of a puzzle. Thanks to the modularity of `ClimaLand` design, many configuration of LSM can be assembled in integrated models. The same functions (`update_aux!`, `exp_tendency!`, etc.) can be used for standalone and integrated models, and an can be stepped  in the same way.

As well as one file:
- ClimaLand.jl: This file is the main Julia module of `ClimaLand.jl` repository. It contains all functions defined in `/src` in a nested way, for example `ClimaLand.X`, `ClimaLand.Soil.X`, 'ClimaLand.Canopy.X`, etc. When a Julia user install and uses ClimaLand via `]add ClimaLand, using ClimaLand`, they are loading those functions, and are ready to use ClimaLand codebase. 
