#include "const.h"

#ifdef DEBUG

#ifdef ansi
#define warn(bool) {if(bool){(void)fprintf(stderr,\
                           "Warning: (%s) true -- file \"%s\", line %d\n",\
                                           # bool ,__FILE__,__LINE__);\
                                         }}

#define fatal(bool) {if(bool){(void)fprintf(stderr,\
                           "fatal error: (%s) true -- file \"%s\", line %d\n",\
                                # bool, __FILE__,__LINE__);\
                                abort(-1) ;\
                              }}
#else
#define warn(bool) {if(bool){(void)fprintf(stderr,\
                           "Warning: (%s) true -- file \"%s\", line %d\n",\
                                           "bool" ,__FILE__,__LINE__);\
                                         }}

#define fatal(bool) {if(bool){(void)fprintf(stderr,\
                           "fatal error: (%s) true -- file \"%s\", line %d\n",\
                                "bool" , __FILE__,__LINE__);\
                                abort(-1) ;\
                              }}
#endif

#else
#define warn(bool)
#define fatal(bool)
#endif

