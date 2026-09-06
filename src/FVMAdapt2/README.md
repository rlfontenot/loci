# FVMAdapt2

FVMAdapt2 is the independently selectable continuation of FVMAdapt. It keeps
the established refinement-plan algorithm while providing an explicit mesh
transition interface for solver-owned data.

The installed components are:

- `fvmadapt2_m.so`: Loci rules for planning and constructing the adapted mesh;
- `libfvmadapt2func`: the supporting tree, geometry, and handoff library; and
- `include/FVMAdapt2`: its public C++ interfaces and refinement-tree types.

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
  depend on private plan encodings or tree keys. State builders and their
  shared declarations remain in `library/mesh_state.h`.

FVMAdapt2 owns old/new topology, geometric overlap, face orientation, node
construction weights, and the communication needed to deliver those
relations. Solvers retain ownership of field semantics, physical bounds,
equation-of-state recovery, boundary policy, and overset connectivity.

## Mesh transition contract

The preferred grid-installation entry point is:

```cpp
setupFVMGridFromContainer(facts, *gridData, cellwts);
```

It installs the mesh and the validated handoff carried by `refinedGridData`.
The legacy overload that accepts the mesh containers separately remains
available, but it installs no transition data.

The cell handoff publishes `refinementDepth`, the number of tree edges from a
current cell to its base cell; `rootCellFileNumber`, the base-VOG cell file
number at the root of that tree; and optional `adaptResult`, which identifies
cells derefined, retained, or refined in the latest transition. A grid
initialized directly from a saved plan has depth and root lineage but no
latest-transition result.

`AMRrefinementMapping` also provides a target-owned `AMRRemapPlan` containing
source and target geometry plus conservative cell contributions. Each
contribution supplies overlap volume and centroid. Its source fraction
distributes an extensive source value; its target fraction assembles a target
average. Components that form a linear bundle may opt into one shared
refinement limiter; existing interpolation calls retain their component-wise
behavior.

For supported transitions, online adaptation carries persistent node, cell,
and face identities in `refinedGridData`, together with node and face remaps.
`setupFVMGridFromContainer` installs these as `nodeId`, `cellId`,
`faceId`, `nodeRemap`, and `faceRemap` facts. The node remap gives
each current node's immediately previous node IDs and geometric weights. The
face remap gives target-owned source-to-target face overlaps with shared
geometry and orientation, and identifies cell-interior faces created or
removed by adaptation. `nodeTransitionReport` and
`faceTransitionReport` state whether the corresponding remap is usable.
Here, source means the mesh before adaptation and target means the resulting
mesh. `FaceOverlap` contains `source` and `target` face IDs, `area`, `centroid`,
and `orientation`. Area is unnormalized; orientation is +1 or -1.
For each target face, `overlaps()` returns its old-face overlaps;
`isCreatedFace()` instead returns the immediately previous cell containing a
new internal face. For an integrated directed flux, the contributor weight is
`overlap.orientation * overlap.area / sourceFaceArea`. A face average uses the target
face area as its denominator. New internal faces have no old-face values to
transfer; their initialization belongs to the caller.

Persistent identities derive from base-mesh entities and refinement-tree paths
rather than generated entity numbers or MPI ownership. They are not Loci
entity numbers or file numbers; `faceId`, `cellId`, and `nodeId` associate
the persistent identities with entities in the installed mesh.

The compact cell plan remains the executable and serialized refinement input.
Transition state retains its accepted leaf paths in a normalized form, so
later transitions do not need to reinterpret topology-specific plan bytes.

On each rank, remap target geometry and rows cover locally owned targets.
Source geometry is the subset needed by the local relation, and value vectors
use the corresponding geometry-array order. Transition-report counts are
reduced over the full MPI job.

Repeated refinement and derefinement handoffs are supported in serial and MPI
for static, planar-faced meshes made of hexahedra, prisms, or a mixture of the
two. If the in-memory state is unavailable, the supported static state can be
replayed from the base VOG and cumulative plan. Plan-start initialization
likewise installs the canonical identities needed by a later transition. The
opaque `MeshState` in `refinedGridData::transitionState` may instead be
carried unchanged into the next in-process adaptation. It retains the accepted
node IDs and coordinates as well as face history, avoiding replay when that
state is available.

Existing `onlineRefineMesh` calls request both remaps. A caller needing only
face transfer can pass `false` as the final `remapNodes` argument to the
overload taking `MeshState`. This skips node-state replay, ancestry records,
and node-remap construction; no `nodeId` fact is installed, the node remap is
empty, and its report is `not_requested`. Node remapping can be enabled on a
later call, which reconstructs the missing source state.
Internally, `collectNodeAncestry` controls recording during coordinate
generation and defaults to false for ordinary mesh-output queries.

General cells can also return face and cell identities, old-face contributors,
and the origins of new internal faces. The tetrahedron regression covers
refinement, retention, and collapse to the root, including MPI empty partitions
and reconstruction of the previous state from its plan. Overlap calculations
require planar convex polygons. A warped internal face can be created, removed,
or retained with an identical vertex loop; its reported area is the magnitude
of its area vector and its center is the wireframe center used by Loci's default
face geometry. Subdivision or merging of warped inherited faces is reported as
unavailable, since no surface-overlap convention has been defined for it.
General-cell node ancestry remains unsupported and is reported independently
of the face handoff.

The public handoff does not define moving-mesh or GCL history, periodic partner
and transform ancestry, overset connectivity, or durable serialization across
process restarts. Those interfaces require metadata from their owning modules.
