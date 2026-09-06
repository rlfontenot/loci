# General-cell face handoff

This test uses the public FVMAdapt2 interface to refine a single tetrahedron,
retain the refined mesh, and coarsen it back to the original cell. It checks
the facts installed in the caller's new grid:

- Boundary faces name their old-face contributors, with positive overlap
  area and the expected orientation. Integrated boundary flux is preserved.
- Each new internal face names the immediately previous cell that contained
  it and has no old-face contributors.
- Retained faces, including warped internal faces, keep their identities and
  receive their full contribution from the same old face.
- Coarsening recombines boundary faces, reports removed internal faces, and
  restores the original face and cell identities.

Both retained history and replay from `currentPlan` are exercised on one, two,
and three MPI ranks. The extra ranks exercise empty cell partitions.
`mpi_*.log` and `debug/` are retained for inspection, including schedules.
The serial log prints each overlap's `source`, `area`, and `orientation`,
plus its source-area fraction, or `sourceCell` for a new internal face.
These IDs correspond to
`faceId` and `cellId`, not the mesh's regenerated entity numbers. The
fraction is overlap area divided by the source face's area; multiply by
`orientation` when transferring a directed quantity.

These tests live on 369 and require a build containing FVMAdapt2. From a
combined checkout with a local install:

```sh
make -C quickTest/FVMAdaptTests/Module/FaceRemap \
  LOCI_BASE="$PWD/loci_install" \
  TEST_BASE="$PWD/quickTest"
```

Clean this directory before changing `LOCI_BASE` to another build.
When running directly from 369, set `LOCI_BASE` to the 390 build or install.
