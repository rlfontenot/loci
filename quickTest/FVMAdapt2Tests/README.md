# FVMAdapt2 tests

All expanded adaptation tests live here on `369-fvmadapt_tests`. The original
upstream `../FVMAdaptTest` suite remains separate and unchanged.

- `Core/Adaptation` checks refinement trees, plan operations, and refinement
  depth, including the edge-ordering and level-refinement regressions.
- The other `Core` groups check conservative transfer and cell, face, and node
  remap contracts.
- `Module` covers XML, parameter files, tags, plan restarts, thin extruded
  meshes, refinement-state facts, and repeated face/node handoffs in MPI.
- `Illustrations` produces optional VTK examples, outside the pass/fail suite.

These tests require the FVMAdapt2 sources and build from 390. From a combined
integration checkout with a local install:

```sh
make -C quickTest FVMAdapt2Tests -j4 LOCI_BASE="$PWD/loci_install"
```

Use `LOCI_BASE="$PWD/OBJ"` to test an uninstalled build. The suite runs 18 test
groups, including serial and two- or three-rank MPI cases. Offline cases use
`marker2`, `refmesh2`, and `refine2`; the original tools still load FVMAdapt.

The checked-in references retain the earlier behavioral baseline: mesh counts,
volume, convexity, and extruded refinement modes. Moving the tests does not
regenerate or relax those references. Face-transition cases also retain logs
and schedules for inspection.
The RefMesh case reuses the upstream suite's mesh and reference files read-only.

The upstream suite is run independently with `make -C quickTest FVMAdaptTest`.
Select the same `LOCI_BASE` explicitly, and clean before changing builds.
