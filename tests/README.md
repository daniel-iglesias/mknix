# MkniX Tests

MkniX uses CMake and CTest. Configure and build from the repository root before running tests:

```sh
cmake -S . -B build
cmake --build build
```

Run the complete suite from the build directory:

```sh
ctest --test-dir build --output-on-failure
```

The legacy tests-directory invocation is also supported:

```sh
ctest --test-dir build/tests --output-on-failure
```

## Unit Tests

The focused C++ unit tests are registered as four CTest tests:

- `mknix_core_unit_tests` exercises classes in `src/core`: Material, Point, Node, shape functions, cells, FEM elements, and Gauss points. `CompBar` is also covered when VTK is enabled. Tests are in [unit/core](unit/core), including [unit/core/material_test.cpp](unit/core/material_test.cpp).
- `mknix_reader_unit_tests` covers top-level reader dispatch and simulation settings, constraint and rigid-body reader settings, flex-body parsing and registration, and nested `SectionReader` fields. The per-class sources are in [unit/reader](unit/reader).
- `mknix_simulation_unit_tests` checks the concrete analysis type identifiers and Simulation static settings in [unit/simulation](unit/simulation).
- `mknix_system_unit_tests` covers System queries, Motion interpolation, force and thermal-load assembly, equal-and-opposite thermal contact, and the distance-constraint penalty response in [unit/system](unit/system).

Build and run only these unit tests with:

```sh
cmake --build build --target mknix_core_unit_tests mknix_reader_unit_tests mknix_simulation_unit_tests mknix_system_unit_tests
ctest --test-dir build -R '^mknix_(core_unit|reader_unit|simulation_unit|system_unit)_tests$' --output-on-failure
```

CTest reports four unit tests passed when all executables succeed. Use `--verbose` to show the individual `[PASS]` cases emitted by the test executables.

## End-to-End Tests

The seven end-to-end tests run `mknixrunner` with input files and compare the generated output against the checked-in reference output. Comparison ignores lines beginning with `GIT_`.

| CTest name | Scenario |
| --- | --- |
| `test1_thermalfluence` | Two-triangle thermal-fluence case |
| `test2_thermalfluence_timefile` | Thermal fluence with a time-dependent input |
| `test3_thermalflux` | Lamella thermal-flux FEM case |
| `test4_thermalflux_RPIM` | Lamella thermal-flux RPIM case |
| `test5_thermalflux_timefile` | Lamella thermal flux with a time-dependent input |
| `test6_thermalflux_2Dmap` | Lamella thermal flux with a 2D map |
| `test7_thermalcontact` | Lamella thermal contact case |

A scenario passes when the generated `.mec` output matches its reference. On failure, CTest displays the captured program output and the comparison error when run with `--output-on-failure`.

## Expected Results

The complete suite contains 11 CTest tests: 4 unit tests and 7 end-to-end tests. A successful full run ends with a summary equivalent to:

```text
100% tests passed, 0 tests failed out of 11
```

The unit-only command reports 4 tests passed. Individual test durations and output formatting may vary by platform and build configuration. The VTK-specific `CompBar` test is compiled only when the project is configured with VTK support.