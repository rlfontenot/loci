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

#ifndef EXT_NAMESPACE
EXT_NAMESPACE=std
#endif

#define HASH_MAP(S,T) EXT_NAMESPACE::hash_map<S,T > 

#endif

#endif
