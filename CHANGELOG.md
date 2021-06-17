# PowerModelsStability.jl Change Log

## staged

- none

## v0.2.1

- Updated README and docs for initial release to Julia Registry
- Fix bug in `add_inverters!` where grounded bus / wye gens were being added even if the data was already kron reduced

## v0.2.0

- Updated to PowerModelsDistribution v0.11.2 (breaking)
- Added `parse_file` function to take both network file and inverter data file paths and parse them appropriately
- Included all keys from inverter data file into the main data structure
- Changed over to new PMD extensions paradigm, and included custom `transform_data_model`, `instantiate_mc_model`, and `solve_mc_opf` functions

## v0.1.0

- Initial release

## v0.0.1
- Initial Commit
- Create basic code structures for future development
- Add the process to generate state space matrix and obtain eigenvalues
