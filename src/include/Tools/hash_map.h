#ifndef LOCAL_HASH_MAP_H
#define LOCAL_HASH_MAP_H

#ifdef USE_MAP_FOR_HASH_MAP

#include <map>
#define HASH_MAP(S,T) std::map<S,T > 

#else

#ifdef EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <hash_map>
#endif


#ifdef GCC_3_1
#define HASH_MAP(S,T) __gnu_cxx::hash_map<S,T > 
#else
#define HASH_MAP(S,T) std::hash_map<S,T > 
#endif

#endif

#endif
