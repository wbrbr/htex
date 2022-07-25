# Htex

This repository provides source code and a demo for the paper ["Htex: Per-Halfedge Texturing for Arbitrary Mesh Topologies"](https://arxiv.org/abs/2207.05618)

## Contents

- `Htex.glsl` is a GLSL library that implements the seamless texture fetch algorithm described in our paper
- The `apps/` directory contains the source code for our Htex viewer `htxviewer`, as well as some utilities:
  - `ptx2htx` converts a Ptex file to Htexthat convert from Ptex to Htex (`ptx2htx`) and from UV mapped meshes to Htex (`uv2htx`)
  - `uv2htx` converts a UV mapped mesh to Htex
  - `aobaker` bakes the ambient occlusion of an input mesh to a Htex file. To build this utility you must toggle the CMake option `BUILD_AO_BAKER`.

You can find the C++ library that we use to read and write Htex files in [this repository](https://github.com/wbrbr/htex-api)