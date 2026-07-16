# BLEND Examples

Examples are organized by purpose:

- `ex01_api` through `ex04_api`: C API examples.
- `ex01_monotone`: command-line `blend monotone` example.
- `ex02_window2d`: command-line `blend window2d` version of `ex02_api`.
- `ex03_window2d`: command-line `blend window2d` South America example.
- `ex04_window3d`: command-line `blend window3d` South America example.
- `ex05_window2d`: command-line `blend window2d` multiple-support example.
- Window-function directories such as `cosine`, `trapezoid`, and `gaussian`:
  one `window1d` example for each function listed by `blend --show-windows`.

Module examples can be run from their own directories:

```sh
cd doc/examples/cosine
./cosine.sh
gnuplot cosine.gp
```

Set `BLEND` to use a specific executable:

```sh
BLEND=/path/to/blend ./cosine.sh
```

The API examples can be built by CMake when `BLEND_BUILD_EXAMPLES` is enabled
in `cmake/ConfigUser.cmake`:

```cmake
set (BLEND_BUILD_EXAMPLES ON)
```

From the repository root:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The API examples can also be run directly from their scripts. If BLEND is
installed, set `BLEND_PREFIX` to the same path used for `CMAKE_INSTALL_PREFIX`:

```sh
cd doc/examples/ex01_api
BLEND_PREFIX=/path/to/blend/install ./ex01_api.sh
```

On macOS, if the dynamic loader cannot find `libblend.dylib`, set:

```sh
export DYLD_LIBRARY_PATH=${BLEND_PREFIX}/lib:${DYLD_LIBRARY_PATH}
```

On Linux, use:

```sh
export LD_LIBRARY_PATH=${BLEND_PREFIX}/lib:${LD_LIBRARY_PATH}
```
