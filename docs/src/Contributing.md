# Contributing

Thank you for contributing to `ClimaLSM`! We encourage
Pull Requests (PRs). Please do not hesitate to ask questions.

## Some useful tips
- When you start working on a new feature branch, make sure you start from
  main by running: `git checkout main`.
- Make sure you add tests
  for your code in `test/` and appropriate documentation in the code and/or
  in `docs/`. All exported functions and structs must be documented.
- When your PR is ready for review, clean up your commit history by squashing
  and make sure your code is current with `ClimateMachine` main by rebasing.

## Continuous integration

After rebasing your branch, you can ask for review. Fill out the template and
provide a clear summary of what your PR does. When a PR is created or
updated, a set of automated tests are run on the PR in our continuous
integration (CI) system.

### Automated testing

Currently a number of checks are run per commit for a given PR.

- `JuliaFormatter` checks if the PR is formatted with `.dev/climaformat.jl`.
- `Documentation` rebuilds the documentation for the PR and checks if the docs
  are consistent and generate valid output.
- `Tests` runs the file `test/runtests.jl`,  using `Pkg.test()`. These are a mix of
  unit tests and fast integration tests.

We use [`bors`](https://bors.tech/) to manage merging PR's in the the `ClimaLSM` repo.
If you're a collaborator and have the necessary permissions, you can type
`bors try` in a comment on a PR to have integration test suite run on that
PR, or `bors r+` to try and merge the code.  Bors ensures that all integration tests
for a given PR always pass before merging into `main`.