# ClimaLand folder structure

[ClimaLand](https://github.com/CliMA/ClimaLand.jl) home directory has 5 main folders:

- docs: contains files to generate [the documentation website](clima.github.io/ClimaLand.jl/).
- experiments: contains simple runs of `ClimaLand` models.
- src: contains the code of `ClimaLand` models.
- ext: contains a package extension for Nueral Snow and an extension with constructors for model paramaters
- lib: contains two packages that extend ClimaLand functionality: ClimaLandSimulations and ClimaLandDashboards
- test: contains [unit tests](https://en.wikipedia.org/wiki/Unit_testing), which are meant to ensure small pieces of ClimaLand source code work as intended before merging pull requests.

and 3 GitHub actions folders. GitHub actions are `.yml` files, which are bash scripts that runs on a remote computer on each git push.

- .buildkite: contains a script building outputs such as figures from experiments and tests folders. These runs are carried out as part of CI and must run without error in order to merge a PR.
- .dev: contains useful tools for developers, such as a format checker for Julia (which is run as part of CI and must pass before a PR can be merged into main).
- .github: contains various scripts, for example, this documentation website is built each time a change is push to `ClimaLand`.

as well as 9 files:

- .gitignore: [commonly used git file](https://git-scm.com/docs/gitignore), contains files, files type, and folders that should be ignored by git.
- LICENSE: License file of `ClimaLand`, you can [read it](https://github.com/CliMA/ClimaLand.jl/blob/main/LICENSE) to learn about legal practice regarding use of `ClimaLand` open source code.
- Project.toml: The Julia programming language requires a `Project.toml` file to create an environment, which specify dependencies of a project as well as its version, name, authors and a unique identifier number (uuid). Every Julia registered package has a `Project.toml` file.
- README.md: This markdown file contains the info that you can read on [ClimaLand GitHub web page](https://github.com/CliMA/ClimaLand.jl)
- logo.svg and logo-white.svg: The logos used in the README
- Artifacts.toml: This file contains `Pkg` artifacts declarations
- NEWS.md: Release notes file with descriptions of new features and changes in each release and main
- NOTICE: This contains the copyright and states what the project is licensed under


## /docs folder

Julia packages are recommended to have a `\docs` folder that builds a standardised documentation following the official documentation generator for Julia: [Documenter.jl](https://documenter.juliadocs.org/stable/).

The folder /docs contains:

- a `/src` folder: It is recommended to put your markdown pages inside this folder. Each markdown file (.md extension text file) is a page accessible through the menu of the documentation. For example, `docs/src/Contributing.md` contains the text you can read on the documentation "contribution guide" menu. The path to this .md file and the name of the menu is set in the `docs/make.jl` file.
- a `make.jl` file: This Julia file contains your documentation website structure. Running this file will build your website pages, you can run it locally, but it is commonly built remotely via .github/workflows/docs.yml to generate the github static page hosted on the [gh-pages branch](https://github.com/CliMA/ClimaLand.jl/tree/gh-pages).

Note: the documentation can have submenu. For example, `APIs` have submenu `ClimaLand` which has many submenu... This structure is built in our current framework via a file `docs/list_of_apis.jl` in that example, which is then included in `docs/make.jl`.

## /experiments folder

The `experiments` folder contains four folders. It contains a folder for `integrated` models and a folder for `standalone` models. These two folders are meant to provide users with simple examples of ClimaLand runs. The files contains meteorological inputs (such as precipitation), values for every parameters, and the domains and timestepper are specified.

For example, `/experiments/integrated/fluxnet/` contains:

- data_tools.jl: Provides utilities for running fluxnet site experiments
- fluxnet_domain.jl: Sets up the domain to run Clima Land on a fluxtower site
- fluxnet_simulation.jl: Contains the site-generic time variables for running `ClimaLand` on
fluxtower sites
- met_drivers_FLUXNET.jl: Construct the drivers
- ozark_pft.jl: Tests running the Ozark site (US-MOz) using plant parameters
defined by plant functional types instead of fully site-specific parameters.
- plot_utils.jl: Plotting utilities for the integrated fluxnet site experiments
- pull_MODIS.jl: provides methods for interacting with the MODIS REST API at a site
- run_fluxnet.jl: Sets up integrated models and runs fluxnet expirement on a given site

and folders for four fluxnet sites. For example `/experiments/integrated/fluxnet/US-MOz` contains:

- US-MOz_parameters.jl: Contains site-specific model parameters for running Clima Land on the Ozark
fluxtower site
- US-MOz_simulation.jl: Contains simulation variables for running Clima Land on the US-MOz
fluxtower site

The `experiments` folder also contains a `benchmarks` folder and a `long_runs` folder. These contain experiments that are run as part of CI on the Caltech cluster.

## /src folder

The `/src` folder contains the source code of `ClimaLand` models. It contains 4 folders:

- shared_utilities: This is a core folder that defines functions and data structures used across all modules and models types of `ClimaLand`. For example, `shared_utilities/models.jl` defines and export the function `make_update_aux` which will be used to create a function which updates the auxiliary parameters, stored in the vector `p`, `shared_utilities/boundary_conditions.jl` defines functions for setting boundary condition for PDE domains, etc.
- standalone: This folder contains standalone models, which are submodels that can be run independently of each other. This is an important aspect of `ClimaLand` code design: to maximize modularity, sub-models can be run alone, and many different methods of the same sub-model can be defined via Julia multiple dispatch. The standalone folder is independent from the integrated folder.
- integrated: This folder contains integrated models. It assembles standalone models together, as one would assemble pieces of a puzzle. Thanks to the modularity of `ClimaLand` design, many configuration of LSM can be assembled in integrated models. The same functions (`update_aux!`, `exp_tendency!`, etc.) can be used for standalone and integrated models, and an can be stepped  in the same way.
- diagnostics: This folder contains default diagnostics methods for different
ClimaLand models. The diagnostics contains metadata such as where diagnostics variables are
stored in those models, what is there names (short, long, standard)
and physical units, with additional comment

As well as two files:

- ClimaLand.jl: This file is the main Julia module of `ClimaLand.jl` repository. It contains all functions defined in `/src` in a nested way, for example `ClimaLand.X`, `ClimaLand.Soil.X`, `ClimaLand.Canopy.X`, etc. When a Julia user install and uses ClimaLand via `]add ClimaLand, using ClimaLand`, they are loading those functions, and are ready to use ClimaLand codebase.
- Artifacts.jl: Contains functions that return the path of the given artifact name in the current context

## /test folder

The `/test` folder contains tests that should