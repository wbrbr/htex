# Htex

This repository provides source code and a demo for the paper ["Htex: Per-Halfedge Texturing for Arbitrary Mesh Topologies"](https://arxiv.org/abs/2207.05618)

## Contents

- `Htex.glsl` is a GLSL library that implements the seamless texture fetch algorithm described in our paper
- The `apps/` directory contains the source code for our Htex viewer `htxviewer`, as well as utilities that convert from Ptex to Htex (`ptx2htx`) and from UV mapped meshes to Htex (`uv2htx`)
- The `lib/` directory contains a C++ library to store and sample Htex files on the CPU. It is a modified version of the Ptex library.