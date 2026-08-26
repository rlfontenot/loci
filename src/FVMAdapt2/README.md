# FVMAdapt2

FVMAdapt2 is the independently selectable continuation of FVMAdapt. It keeps
the established refinement-plan algorithm while providing an explicit mesh
transition interface for solver-owned data.

The installed components are:

- `fvmadapt2_m.so`: Loci rules for planning and constructing the adapted mesh;
- `libfvmadapt2func`: the supporting tree, geometry, and handoff library; and
- `include/FVMAdapt2`: its public and internal C++ interfaces.

FVMAdapt and FVMAdapt2 intentionally provide many of the same Loci facts and
C++ entry points. A program must select and load one implementation, never
both. Keeping separate module, library, and include names lets existing
applications continue using FVMAdapt while applications that need the richer
transition contract opt into FVMAdapt2.

## Source organization

- The top-level `.loci` files expose scheduler-visible plans, relations,
  constraints, numbering, and mesh construction stages.
- `library/` contains local refinement-tree algorithms and output adapters.
- Public transition types live under `include/FVMAdapt2`; solvers should not
  depend on private plan encodings or tree keys.

FVMAdapt2 owns old/new topology, geometric overlap, face orientation, node
construction weights, and the communication needed to deliver those
relations. Solvers retain ownership of field semantics, physical bounds,
equation-of-state recovery, boundary policy, and overset connectivity.

The public face and node remap classes currently define and validate the
intended transition contract. The online adaptation path does not yet populate
those relations; consumers must not treat their presence in the headers as a
live ancestry handoff.
