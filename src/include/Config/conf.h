#if 0
#if __GNUC__ == 3 
# if __GNUC_MINOR__ == 2
#  define GCC_3_2
# else 
#  define GCC_3_0
# endif
#endif
#if __GNUC__ == 2
# if __GNUC_MINOR__ == 95
#  define GXX_FIXES
# endif
#endif

#ifdef GCC_3_0
#define EXT_HASH_MAP
#endif

#ifdef GCC_3_2
#define EXT_HASH_MAP
#define EXT_NAMESPACE __gnu_cxx
#define NO_OFFSETOF
#endif
#endif

