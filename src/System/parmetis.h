#define PARMETIS_MAJOR_VERSION        3
#define PARMETIS_MINOR_VERSION        0

typedef int idxtype;


/*****************************/
/* NEW Partitioning Routines */
/*****************************/
/* kmetis.c */
void ParMETIS_V3_PartKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, float *, float *, int *, int *, idxtype *, MPI_Comm *);

/* mmetis.c */
void ParMETIS_V3_PartMeshKway(idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int *, idxtype *, MPI_Comm *);

/* gkmetis.c */
void ParMETIS_V3_PartGeomKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, float *, float *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_V3_PartGeom(idxtype *, int *, float *, idxtype *, MPI_Comm *);


/************************/
/* Adaptive subroutines */
/************************/

/* ametis.c */
void ParMETIS_V3_AdaptiveRepart(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, float *, float *, float *, int *, int *, idxtype *, MPI_Comm *);

/* rmetis.c */
void ParMETIS_V3_RefineKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, float *, float *, int *, int *, idxtype *, MPI_Comm *);


/****************************/
/* Mesh to Dual subroutines */
/****************************/
/* mesh.c */
void ParMETIS_V3_Mesh2Dual(idxtype *, idxtype *, int *, int *, int *, idxtype **, idxtype **, MPI_Comm *);


/************************/
/* Ordering subroutines */
/************************/
/* ometis.c */
void ParMETIS_V3_NodeND(idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, MPI_Comm *);

/*************/
/* pspases.c */
/*************/
void ParMETIS_SerialNodeND(idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, MPI_Comm *);


/***************************************/
/* backwards compatibility subroutines */
/***************************************/
/* kmetis.c */
void ParMETIS_PartKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void PARKMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
/* rmetis.c */
void ParMETIS_RefineKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void PARRMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
/* diffuse.c */
void ParMETIS_RepartLDiffusion(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_RepartGDiffusion(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void PARUAMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PARDAMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
/* scremap.c */
void ParMETIS_RepartRemap(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_RepartMLRemap(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
/* gmetis.c */
void ParMETIS_PartGeomKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_PartGeom(idxtype *, int *, float *, idxtype *, MPI_Comm *);
void ParMETIS_PartGeomRefine(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *, MPI_Comm *);
void PARGKMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);
void PARGRMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);
void PARGMETIS(idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);
/* ometis.c */
void ParMETIS_NodeND(idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, MPI_Comm *);
void PAROMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);


