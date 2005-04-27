
#if defined(__GNUC__)

#if defined(__INTEL_COMPILER)
// Intel Compiler

#else

// GCC Compilers
#if (__GNUC__==3) && (__GNUC_MINOR__==0)
#define EXT_HASH_MAP
#endif

#if (__GNUC__==3) && (__GNUC_MINOR__!=0)
#define EXT_HASH_MAP
#define EXT_NAMESPACE __gnu_cxx
#define NO_OFFSETOF
#endif

#if (__GNUC__>=4)
#define EXT_HASH_MAP
#define EXT_NAMESPACE __gnu_cxx
#define NO_OFFSETOF
#endif

#endif
#else

#if defined(__sgi)
// SGI Compiler
#define NO_CSTDLIB
#define NO_CMATH
//#define MPI_NO_CPPBIND
#endif

#endif


#if defined(LOCI_SYS_IRIX64)
#undef SGI
#define SGI
#endif

#if defined(LOCI_SYS_Linux)
#undef LINUX
#define LINUX
#endif

#if defined(LOCI_SYS_SunOS)
#undef SPARC
#define SUN
#define SPARC
#endif

