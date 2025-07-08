# GreenX Library - TimeFrequency

This library provides optimal quadrature grid points and weights for imaginary time-frequency transforms, commonly 
occurring in MP2, RPA and Green's function methods. Optimisation is performed with the Minimax procedure, minimising the
maximum error of the quadrature. This typically results in an error that is more evenly distributed across the interval 
of interest. Grids are provided with n-points ranging from 6 to 34, and for different values of the transition energy 
ratios Rm (on average 15 R-values for each grid point).

For additional details, please refer to the corresponding JOSS paper, included [here](../JOSS).

## Building

With CMake, change to the GreenX root, then type:

```bash
mkdir build && cd build
cmake ../
make -j 
make install 
```

## Running the Tests

Application tests are run with the pytest framework. Having installed `pygreenx`
(see the top-level [README](../README.md)) change to `<GX_ROOT>/<BUILD_DIR>` and type `ctest`. 
Additionally, one can change to the test folder and explicitly run the pytest 
command from there:

```bash
cd <GX_ROOT>/<BUILD_DIR>/test/time-frequency
pytest -s 
# or
pytest -s test_name.py
```

## Calling the MiniMax Library

Information on how to use the library is always up-to-date on the [GreenX website](nomad-coe.github.io/greenX/).

## For Developers

### Adding a New Test

Write a `test_name.py` and `test_name.f90` in `test/`. 
`test_name.f90` should accept some input parameter/s, make a library call and
write the result to file. `test_name.py` should provide the input parameters
(the simplest way is to write to file), execute the binary built from 
`test_name.f90`, parse the result and make the assertions. This implies the
reference results should also be stored in `test_name.py`.

`test_name.py` *must* include the fixture:

```python
@pytest.fixture()
@pytest.mark.usefixtures("get_binary", "greenx_build_root")
def fortran_binary(get_binary, greenx_build_root):
    name = 'gx_tabulate_grids.exe'
    _binary = get_binary(name)
    assert _binary is not None, f'{name} cannot be found in {greenx_build_root}'
    print(f'Binary source: {_binary}')
    return _binary
```

where the binary `name` should be consistent with whatever the developer has
called the binary in the CMakeLists.txt

One notes that the working directory for all tests is `tmp_path`, which is 
switched to using a fixture defined in `GX-TimeFrequency/test/conftest.py`, 
for the duration of the test execution. The `tmp_path` fixture will provide a 
temporary directory unique to the (py)test invocation, created in the base 
temporary directory.

A new test can be added to this `CMakeLists.txt` in the following way:

```cmake
set(LIBS_FOR_TESTING)
list(APPEND LIBS_FOR_TESTING LibGXMiniMax)

set(TEST_NAME "test_gx_minimax_grid")
add_app_test(TEST_NAME TEST_TARGET_DIR LIBS_FOR_TESTING)
```

Calling `add_app_test` should ensure that the directory level structure in
the build folder is consistent with example python driver. For multiple tests, 
one can place the call to `add_app_test` in a loop, and supply different test names.

All the tests can be run by switching to the `build/` directory and typing `ctest`,
however they can also be run in the `test/` directory of the source by typing
`pytest -s` as long as the 'ENV VAR' to the build directory has been defined:

```bash
export GX_BUILD_DIR=<path/2/build/>
```
