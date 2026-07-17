# fastwindingnumber

Native C++17 implementation of fast winding number queries for oriented triangle meshes.
The project uses only the C++ standard library.

## Build

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## Command Line

```bash
./build/fwn_query data/cube.obj data/sample_points.xyz --fast
```

Input mesh is a triangular or polygonal OBJ file. Polygonal faces are triangulated with a fan.
Point files contain one query point per line:

```text
x y z
```

Output format:

```text
x y z winding_number
```

With `--classify`, a fifth column is added. It is `1` when `abs(winding_number) >= threshold`,
otherwise `0`.

## Benchmark

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
./build-release/fwn_benchmark --queries 4096 --repeats 5 --beta 5.0 --expansion-order 2
```

The benchmark builds a procedural torus mesh, samples random query points inside its bounding box,
and reports exact-vs-fast timing plus worst and mean absolute error.

## Notes

- Correct inside/outside classification requires consistently oriented input surfaces.
- `--exact` evaluates all triangle solid angles directly.
- `--fast` uses a BVH. Far nodes are approximated by an aggregated area-vector multipole expansion,
  while near nodes fall back to exact triangle solid angles.
- `beta` controls the speed/accuracy tradeoff. Larger values visit more nodes and are more accurate.
  The default is `3.5`, which favors accuracy over the original coarse setting.
- `--expansion-order 2` is the default and most accurate approximation mode. `0` reproduces the
  coarse single-dipole approximation, and `1` is a middle ground.
