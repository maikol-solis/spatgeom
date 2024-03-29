# Changelog

All notable changes to this project will be documented in this file.

## [0.3.0] - 2023-04-26

### Bug Fixes

- Remove old code of geometric sensitivity
- Change data_frame_triangles to geom_indices
- Slice the triangles df instead of subset it.
- Correct print and plot functions with the new data structure.
- Return triangle_list in estimate_envelope.
- Use ribbons instead of lines for the envelope.
- Plot envelope ribbon min and max without NAs.
- Force to finite values in envelope ribbon.
- Documentation and import checks.
- Tweak plot styles (in derivative).

### Documentation

- Add package badges.

### Features

- Add linear and donut examples in a separate `data.R` file.

### Miscellaneous Tasks

- Add pre-commit config

### Refactor

- Move envelope loop to its own function

## [0.2.0] - 2023-02-14

### Bug Fixes

- Change alphastats to spatgeom in manual

### Features

- Include envelope plots

## [0.1.0] - 2022-07-02

This is the initial release for `spatgeom`. 

### Bug Fixes

- Fix parameter names to make the package more modular.
- Fix print function to a S3 method.  
- Fix github urls repository.
- Fix miscellaneous things for CRAN checks. 

### Documentation

- Add Roxygen documentation on all relevant functions.
- Add simple enough examples under the `@examples` tag to pass all the CRAN
  tests.

### Features

- Add new parameter for `mc_cores`. Now, the main function could run the
  estimation in parallel cores.

### Refactor

- Delete old barcode and plots implementations.
- Remove unused code overall all the package.  
- Split plot functions into `plot_curves` and `plot_alphashapes`.
