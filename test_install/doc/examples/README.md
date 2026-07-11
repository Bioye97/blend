# BLEND Examples

Examples are organized in numbered directories:

```text
ex01, ex02, ...
```

Each example directory contains:

- `exNN.c`: C source code for the example
- `exNN.gp`: optional gnuplot script

When an example is run, it writes an `exNN.txt` data file in the current
directory. The gnuplot script writes `exNN.png`. The 3D example uses a
color-coded point cloud to show the taper weights in volume.

Current examples:

- `ex01`: symmetric cosine taper window of length `n = 100`
- `ex02`: symmetric trapezoid taper window of length `n = 100`
- `ex03`: 2D cosine taper over a triangular support on a 100 by 100 grid
- `ex04`: 2D cosine taper over a square support on a 100 by 100 grid
- `ex05`: 2D cosine taper over a pentagon support on a 100 by 100 grid
- `ex06`: 2D cosine taper over a hexagon support on a 100 by 100 grid
- `ex07`: 2D cosine taper over a heptagon support on a 100 by 100 grid
- `ex08`: 2D cosine taper over an octagon support on a 100 by 100 grid
- `ex09`: 2D cosine taper over a nonagon support on a 100 by 100 grid
- `ex10`: 2D cosine taper over a decagon support on a 100 by 100 grid
- `ex11`: 2D cosine taper over a 4-pointed isotoxal-star support on a 100 by
  100 grid
- `ex12`: 3D cosine taper over the same 4-pointed isotoxal-star support as
  `ex11`, using a 100 by 100 by 25 grid

The 1D examples write:

```text
index weight
```

The 2D examples write:

```text
x y weight
```

The 3D examples write:

```text
x y z weight
```

## Build Examples With CMake

Enable examples in `cmake/ConfigUser.cmake`:

```cmake
set (BLEND_BUILD_EXAMPLES ON)
```

From the repository root, build BLEND:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

To run an example from its example directory:

```sh
cd ../doc/examples/ex01
../../../build/ex01
```

Replace `ex01` with another example name, such as `ex02`.

## Compile Directly From the Source Tree

From the repository root:

```sh
EX=ex01
BLEND_SOURCES="src/blend_boundary.c src/blend_contribution.c src/blend_polygon.c src/blend_window.c"
cc ./doc/examples/${EX}/${EX}.c ${BLEND_SOURCES} -I src -o ./doc/examples/${EX}/${EX} -lm
```

Then run it:

```sh
cd ./doc/examples/${EX}
./${EX}
```

## Compile With an Installed Shared Library

If BLEND is installed, set `BLEND_PREFIX` to the same path used for
`CMAKE_INSTALL_PREFIX` in `cmake/ConfigUser.cmake`:

```sh
export BLEND_PREFIX=/path/to/blend/install
```

Then compile from the repository root:

```sh
EX=ex01
cc ./doc/examples/${EX}/${EX}.c \
  -I ${BLEND_PREFIX}/include/blend \
  -L ${BLEND_PREFIX}/lib \
  -lblend \
  -o ./doc/examples/${EX}/${EX}
```

On macOS, if the dynamic loader cannot find `libblend.dylib`, set:

```sh
export DYLD_LIBRARY_PATH=${BLEND_PREFIX}/lib:${DYLD_LIBRARY_PATH}
```

On Linux, use:

```sh
export LD_LIBRARY_PATH=${BLEND_PREFIX}/lib:${LD_LIBRARY_PATH}
```

## Plot With Gnuplot

From an example directory after running the example:

```sh
gnuplot ex01.gp
```

Replace `ex01.gp` with the plotting script for the example you ran.
