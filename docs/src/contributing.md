# Contributing

Thank you for contributing to `ClimaLand`! We encourage Pull Requests (PRs).
Please do not hesitate to ask questions, or to open issues if something seems amiss
or you'd like a new feature.

## Some useful tips
- When developing code it's best to work on a branch off of the most recent main.
This can be done by running the following commands, where "initials" corresponds to the first and last initial of the person starting the branch.
```
git checkout main
git pull
git checkout -b initials/branch_name
```

- Make sure you add tests for your code in `test/`, appropriate documentation in `docs/`,
  and descriptive inline comments throughout the code.
  All exported functions and structs must have docstrings.
- When your PR is ready for review, clean up your commit history by squashing to 1 commit per PR
  and make sure your code is current with `ClimaLand.jl` main by rebasing.

## Continuous integration

After rebasing your branch, you can ask for review. Fill out the template and
provide a clear summary of what your PR does. When a PR is created or
updated, a set of automated tests are run on the PR in our continuous
integration (CI) system.

ClimaLand.jl's continuous integration contains multiple automated tests,
which are described below. All of these must pass for a PR to be eligible
to merge into the main branch.

### Unit testing

Unit tests are defined in the `test/` folder, and are all listed in the file
`test/runtests.jl`, which allows them to easily be called together.
In our CI, these tests are run via Github Actions on a variety of operating
systems and Julia versions. We also use downgrade tests to check the lower limits
of our Julia package compat bounds, and downstream tests to verify compatibility
with ClimaCoupler.jl.

### Buildkite pipeline

The buildkite pipeline contains a variety of ClimaLand simulations,
which span the complexity of our models and domains.
Some of these simulations have explicit checks, for example comparing
to empirical solutions or output from external codebases, while others
test that the simulation can run without error.

This pipeline also runs our unit test suite on GPU and with MPI,
to ensure all of our source code is compatible with those setups.

### Formatting check

The `JuliaFormatter` test checks if the PR is correctly formatted according
to `.dev/climaformat.jl`. To ensure that this is the case, we recommend running
`julia --project=.dev .dev/climaformat.jl` on each PR before requresting review.
This command will apply the formatter and any changes it detects are necessary.

### Documentation

The `Documentation` test rebuilds the documentation for the PR and checks if the docs
are consistent and generate valid output.
