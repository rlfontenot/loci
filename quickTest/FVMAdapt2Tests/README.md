# FVMAdapt2 transition tests

These tests belong to `369-fvmadapt_tests`. They require the FVMAdapt2 sources
and build from `390-fvmadapt-rule-syntax-update`; run them in the combined
integration checkout, not against the older implementation on 369 alone.

After building and installing the combined checkout, run from its root:

```sh
make -C quickTest/FVMAdapt2Tests -j4 \
  LOCI_BASE="$PWD/loci_install" TEST_BASE="$PWD/quickTest"
```

The suite checks cell, face, and node remaps, repeated refinement and
derefinement, plan replay, and the installed refinement-state facts. Module
tests include serial and two- or three-rank MPI runs. `TestResults` summarizes
the nine test groups; face-transition cases retain logs and schedules under
their case directories for inspection.

`../FVMAdaptTests/Module/FaceRemap` separately exercises the public general-cell
face handoff using only the installed FVMAdapt2 interface.
