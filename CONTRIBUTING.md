# Developer Guide

In this guide we aim to communicate the various code standards expected for this package.

## Documentation

Documentation should be included for all new publically exported additions to the code base

### Examples

In the case of new functionality, like helper functions for data manipulation or model debugging, or changes to building or running models, notebooks should be utilized to demonstrate these features, and are required to exist in the `/examples` directory.

## Style Conventions

In general, the following conventions should be adhered to when making changes or additions to the code base. These conventions should include any conventions applied across the InfrastructureModels ecosystem specific to power engineering (_i.e_ conventions from InfrastructureModels, PowerModels, PowerModelsRestoration, etc.).

### Functions

Function additions should meeting the following criteria:

- All functions should be clearly named, without abbreviations, and with underscores between words, _e.g._ `parse_file` or `constraint_bus_voltage_magnitude`; in Python this is known as [`lower_case_with_underscores`](https://legacy.python.org/dev/peps/pep-0008/#descriptive-naming-styles). The exception to the abbreviate rule is cases where abbreviations would be expected in the modeling of power systems.
- All functions that are not prepended by an underscore `_` will be exported by default (_i.e._ when a user uses `using`). Public functions should have a detailed docstring instructing on usage
- All functions that modify data in place should end with an exclamation point `!` and the function input that is being modified should be the first argument (or first arguments in the case where multiple inputs are being modified in place). The exceptions to this rule are constraint and variable creation functions (_i.e._ those functions related to JuMP model creation), which do not include the exclaimation point
- All function arguments, including keyword arguments, should have their types specified.
- Private functions, _i.e._ those intended to be for internal use only, should follow the same descriptive naming conventions as functions exported by default, and should always include docstrings to describe their purpose.
- Functions should be separated by two blank lines

```julia
"this function demonstrates how an internal, in-place data altering function should be defined"
function _concise_descriptive_name!(data::Dict{String,<:Any}, a::Real, b::Vector{<:Real}, c::Matrix{<:Complex}; d::Bool=false, e::Vector{<:Function}=Vector{Function}([]))
end
```

### Types & Enums

When specifying types, _i.e._ when specifying the type of a function argument, or creating enums, these guidelines are recommended:

- Prefer to use `Vector{T}` instead of `Array{T,1}`
- Prefer to use `Matrix{T}` instead of `Array{T,2}`

### Constants

Whenever possible, `const` should be used to eliminate unnecesary re-evaluations of code, and every `const` should have a docstring, whether internal or public.

Currently, all phase-aware functions use `mc`, but this is subject to change in the future as we refactor. If the function is not multiphase specific, these are not needed in the function name.

### Metaprogramming

In general, it is better to avoid metaprogramming patterns, like creating functions algorithmically, in order to aid in the debugging of code. Metaprogramming can create significant challenges in interpreting stacktraces upon errors.

### Markdown

Markdown files should be properly formatted, particularly when including tables. Developers are encouraged to use [markdownlint](https://github.com/markdownlint/markdownlint) and a markdown formatter (such as in VSCode).

## File Structure

It is important that new functions, variables, constraints, etc. all go into appropriate places in the code base so that future maintenance and debugging is easier. Pay attention to the current file structure and attempt to conform as best as possible to it. In general

- `/src/core` contains the core logic of the package, including variable creation and constraint templates, _i.e._ things that are agnostic to the formulation
- `src/form` contains formulation specific variable and constraint functions, organized under separate files for different formulations
- `src/io` contains all of the tools to parse and save files, in particular all of the logic necessary to parse dss files and output json files
- `src/prob` contains all problem specifications
- `docs/src` contains all source markdown files for the documentation

## Dependencies (Project.toml)

All new dependencies should be carefully considered before being added. It is important to keep the number of external dependencies low to avoid reliance on features that may not be maintained in the future. If possible, Julia Standard Library should be used, particularly in the case where reproducing the desired feature is trivial. There will be cases where it is not simple to duplicate a feature and subsequently maintain it within the package, so adding a dependency would be appropriate in such cases.

All new dependencies are are ultimately approved should also include an entry under `[compat]` indicating the acceptable versions (Julia automerge requirement). This includes test-only dependencies that appear under `[extras]`

The `Manifest.toml` __should not__ be included in the repo.

## Pull Requests

All pull requests should be reviewed by a core developer, and may include a review by a subject matter expert if the area of the PR is outside that of one of the core developers. In that case, the core developers will primarily review style and design, rather than substance.

Every PR should strive to meet the following guidelines.

### PR Title

- Should be concise and clear, describing in a phrase the content of the PR
- Should include a prefix that describes the primary type of the PR
  - ADD: feature addition
  - FIX: bugfix
  - REF: refactor
  - UPD: updates to code for _e.g._ version bumps of dependencies
  - STY: style changes, no changes to function names, added features, etc.
  - DOC: documentation-only additions/changes
  - RM: dead code removal

### PR Body

- If the change is breaking, it should be clearly stated up front
- The purpose of this PR should be clearly stated right away
- Major changes / additions to the code should be summarized. In the case where a refactor was performed, the name changes of public functions should be documented in the body of the PR
- Any associated Issues should be referenced in the body of the PR, and it is accepted/encouraged to use Closes #XX to automatically close Issues after the PR is merged

### PR Code

- An entry should be added to CHANGELOG.md for every PR
- Documentation should be updated (See [Documentation](##Documentation) section above for guidelines)
- Unit tests should be added. In the case where existing unit tests were altered, an explanation for the change must be included
- Code should be rebased to the latest version of whatever branch the PR is aimed at (no merge conflicts!)

## Versions

We follow the Semantic Versioning ([SemVer](https://semver.org/)) convention of `Major.minor.patch`, where `Major` indicates breaking changes, `minor` indicates non-breaking feature additions, and `patch` indicates non-breaking bugfixes.

Currently, because `Major==0`, `minor` indicates breaking changes and `patch` indicates any non-breaking change, including both feature additions and bugfixes. Once this package reaches `v1.0.0`, we will adhere strictly to the SemVer convention.

## Branch Management

The `master` branch is a [protected](https://help.github.com/en/github/administering-a-repository/about-protected-branches) branch, meaning that its history will always be contiguous and can never be overwritten.

Release candidate branches of the format `vM.m.0-rc` are also protected branches. These branches will contain only breaking changes and will not be merged into master until a new version is ready to be tagged. Pull requests including breaking changes should be directed into the next release candidate branch available, _e.g._ if the current version of the package is `v0.9.0`, the next release candidate branch will be `v0.10.0-rc`.

Pull requests that include only non-breaking changes can be merged directly into `master` once approved, and in the case of merge conflicts arising for release candidate branches, the `-rc` branch will need to be updated to include the latest `master`.

Pull requests will generally be merged using [squash and merge](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-request-merges#squash-and-merge-your-pull-request-commits) into the branch they are aimed at, with the exception of release candidate branches, which generally be merged using [rebase and merge](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-request-merges#rebase-and-merge-your-pull-request-commits) into master.
