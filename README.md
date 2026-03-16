# C++ Recursive Ray Tracer

![C++17](https://img.shields.io/badge/C++-17-blue.svg)
![Rendering](https://img.shields.io/badge/Rendering-Offline-orange)
![License](https://img.shields.io/badge/license-MIT-green)

> A physically based, recursive ray tracer built entirely from scratch in C++ — no external graphics APIs.

---

## Rendered Output

Each image below was produced by the engine from a plain-text scene description file.

<p align="center">
  <img src="images/test1.png" width="45%" title="Sphere scene with Blinn-Phong shading">
  <img src="images/test2.png" width="45%" title="Reflective and refractive materials">
</p>
<p align="center">
  <img src="images/test3.png" width="45%" title="Mixed geometry with shadows">
  <img src="images/earth.png" width="45%" title="Spherical texture mapping — Earth">
</p>

---

## Features

| Feature | Description |
|---|---|
| Recursive ray tracing | Up to 10 bounce depth for reflection and refraction |
| Blinn-Phong illumination | Ambient, diffuse, and specular components per light |
| Fresnel reflectance | Schlick approximation for angle-dependent reflectivity |
| Snell's Law refraction | Physically accurate transmission with total internal reflection |
| Shadow computation | Shadow rays against all occluders with transparency support |
| Geometry | Ray–sphere and ray–triangle intersections |
| Triangle meshes | Flat and smooth (interpolated normal) shading |
| Texture mapping | Spherical UV for spheres, barycentric interpolation for triangles |
| Scene format | Custom `.txt` scene description files |
| Output | Uncompressed PPM image files |

---

## Project Structure

```
Ray-Tracing/
├── 01-ray-casting/           # Basic ray casting — camera, rays, FOV
├── 02-illumination-shadows/  # Blinn-Phong shading + shadow rays
├── 03-mesh-textures/         # Triangle meshes, smooth normals, texture mapping
├── 04-reflection-refraction/ # Full ray tracer — Fresnel, Snell's Law, recursion
├── modern/                   # Modernized C++17 rewrite (Vec3 type, constexpr, std::filesystem)
└── images/                   # Reference renders
```

Each numbered stage is self-contained and builds progressively on the previous one.

---

## Build & Run

### Requirements
- GCC 7+ or Clang 6+ with C++17 support
- GNU Make

### Original implementation (04-reflection-refraction)
```bash
cd 04-reflection-refraction
make
./out test1.txt
```

### Modernized version (modern/)
```bash
cd modern
make
./ray_tracer ../04-reflection-refraction/test1.txt
```

Both executables take a `.txt` scene file and write a `.ppm` image to the same directory.

---

## Scene File Format

Scene files define the camera, lights, geometry, and materials in plain text:

```
eye       0 2 4
viewdir   0 -0.7 -1
updir     0 1 0
vfov      90
imsize    800 800
bkgcolor  0.2 0.2 0.2

light     0 -1 -1  0  1 1 1
mtlcolor  0.1 0.3 0.8  1 1 1  0.1 0.7 0.5 32  1.0 1.0

sphere    0 0 -3  1
```

| Keyword | Description |
|---|---|
| `eye` | Camera position |
| `viewdir` | Viewing direction |
| `vfov` | Vertical field of view (degrees) |
| `imsize` | Output image width and height |
| `bkgcolor` | Background color (RGB 0–1) |
| `light` | Position, type (0=directional, 1=point), RGB color |
| `mtlcolor` | Diffuse RGB, specular RGB, ka kd ks shininess, alpha, eta |
| `sphere` | Center XYZ, radius |
| `v / vn / vt` | Vertex, normal, UV coordinate |
| `f` | Triangle face (v/t/n indices, 1-based) |
| `texture` | Path to a PPM texture file |

---

## Modernized Rewrite

The `modern/` folder contains a C++17 rewrite of the final renderer with the following improvements over the 2019 original:

- **`Vec3` math type** with operator overloads, replacing raw C-style float arrays and manual `Cross`/`Dot`/`Normalize` calls throughout
- **`constexpr` constants** replacing `#define PI` and `#define epsilon`
- **`std::optional<HitRecord>`** for intersection results, eliminating sentinel values like `100000`
- **`std::filesystem`** for output path construction
- **`std::clamp`** for pixel value clamping
- **Named, purpose-clear structs** — `Ray`, `Material`, `Sphere`, `Light`, `Face`, `Texture`, `HitRecord`
- **Möller–Trumbore** ray–triangle intersection replacing the plane-then-barycentric approach
- **Removed all commented-out debug code**
- Compiles cleanly under `-std=c++17 -Wall -Wextra`
