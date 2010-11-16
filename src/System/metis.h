//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef METIS_H
#define METIS_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#define IDXTYPE_INT
#ifdef IDXTYPE_INT
typedef int idxtype ;
#else
typedef short idxtype ;
#endif

#ifdef __cplusplus
extern "C" {
#endif
/* estmem.c */
void METIS_EstimateMemory(int *, idxtype *, idxtype *, int *, int *, int *);

void METIS_PARTGRAPHRECURSIVE(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphrecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphrecursive_(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphrecursive__(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void METIS_WPARTGRAPHRECURSIVE(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphrecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphrecursive_(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphrecursive__(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void METIS_PARTGRAPHKWAY(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphkway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphkway_(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphkway__(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void METIS_WPARTGRAPHKWAY(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphkway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphkway_(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphkway__(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void METIS_EDGEND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_edgend(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_edgend_(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_edgend__(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void METIS_NODEND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodend(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodend_(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodend__(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void METIS_NODEWND(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodewnd(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodewnd_(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodewnd__(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void METIS_PARTMESHNODAL(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void metis_partmeshnodal(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void metis_partmeshnodal_(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void metis_partmeshnodal__(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void METIS_PARTMESHDUAL(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void metis_partmeshdual(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void metis_partmeshdual_(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void metis_partmeshdual__(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void METIS_MESHTONODAL(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtonodal(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtonodal_(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtonodal__(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void METIS_MESHTODUAL(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtodual(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtodual_(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtodual__(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void METIS_ESTIMATEMEMORY(int *, idxtype *, idxtype *, int *, int *, int *);
void metis_estimatememory(int *, idxtype *, idxtype *, int *, int *, int *);
void metis_estimatememory_(int *, idxtype *, idxtype *, int *, int *, int *);
void metis_estimatememory__(int *, idxtype *, idxtype *, int *, int *, int *);
void METIS_MCPARTGRAPHRECURSIVE(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);
void metis_mcpartgraphrecursive(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);
void metis_mcpartgraphrecursive_(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);
void metis_mcpartgraphrecursive__(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);
void METIS_MCPARTGRAPHKWAY(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);
void metis_mcpartgraphkway(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);
void metis_mcpartgraphkway_(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);
void metis_mcpartgraphkway__(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);
void METIS_PARTGRAPHVKWAY(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);
void metis_partgraphvkway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);
void metis_partgraphvkway_(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);
void metis_partgraphvkway__(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);
void METIS_WPARTGRAPHVKWAY(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);
void metis_wpartgraphvkway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);
void metis_wpartgraphvkway_(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);
void metis_wpartgraphvkway__(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);

/* kmetis.c */
void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void METIS_WPartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 

/* kvmetis.c */
void METIS_PartGraphVKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);
void METIS_WPartGraphVKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);

/* mesh.c */
void METIS_MeshToDual(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void METIS_MeshToNodal(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void GENDUALMETIS(int, int, int, idxtype *, idxtype *, idxtype *adjncy);
void TRINODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);
void TETNODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);
void HEXNODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);
void QUADNODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);

/* meshpart.c */
void METIS_PartMeshNodal(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void METIS_PartMeshDual(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);

/* mkmetis.c */
void METIS_mCPartGraphKway(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);


/* mpmetis.c */
void METIS_mCPartGraphRecursive(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);
void METIS_mCHPartGraphRecursive(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);
void METIS_mCPartGraphRecursiveInternal(int *, int *, idxtype *, idxtype *, float *, idxtype *, int *, int *, int *, idxtype *);
void METIS_mCHPartGraphRecursiveInternal(int *, int *, idxtype *, idxtype *, float *, idxtype *, int *, float *, int *, int *, idxtype *);
/* ometis.c */
void METIS_EdgeND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void METIS_NodeND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void METIS_NodeWND(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 

/* parmetis.c */
void METIS_PartGraphKway2(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void METIS_WPartGraphKway2(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void METIS_NodeNDP(int, idxtype *, idxtype *, int, int *, idxtype *, idxtype *, idxtype *);
void METIS_NodeComputeSeparator(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *); 
void METIS_EdgeComputeSeparator(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *); 

/* pmetis.c */
void METIS_PartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void METIS_WPartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
#ifdef __cplusplus
}
#endif
#endif
