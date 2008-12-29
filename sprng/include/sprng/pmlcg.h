
#ifndef _pmlcg_h
#define _pmlcg_h

#ifndef ANSI_ARGS
#ifdef __STDC__
#define ANSI_ARGS(args) args
#else
#define ANSI_ARGS(args) ()
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

int pmlcg_get_rn_int ANSI_ARGS((int *igenptr));
float pmlcg_get_rn_flt ANSI_ARGS((int *igenptr));
double pmlcg_get_rn_dbl ANSI_ARGS((int *igenptr));
int *pmlcg_init_rng ANSI_ARGS((int rng_type,  int gennum, int total_gen,  int seed,
			  int mult));
int pmlcg_spawn_rng ANSI_ARGS((int *igenptr, int nspawned, int ***newgens, int checkid) );
int pmlcg_get_seed_rng ANSI_ARGS((int *genptr));
int pmlcg_free_rng ANSI_ARGS((int *genptr));
int pmlcg_pack_rng ANSI_ARGS(( int *genptr, char **buffer));
int *pmlcg_unpack_rng ANSI_ARGS(( char *packed));
int pmlcg_print_rng ANSI_ARGS(( int *igen));


#ifdef __cplusplus
}
#endif


#endif
