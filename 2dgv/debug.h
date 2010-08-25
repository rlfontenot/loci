#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#include <stdlib.h>

#ifdef DEBUG

/* 
   When DEBUG is defined, warn(bool) and fatal(bool) print warning messages 
   if the boolean argument they are passed evaluates to a true value.  
   fatal(bool) terminates the program after printing the message.  
   If DEBUG is not defined, these macros do not insert any code.
   
   These routines are used to provide some level of assurance about program
   correctness and should not be interpreted as a means of providing warnings
   or error messages to the user.
   */

#define DEBUG_TYPE_NAME(v) (typeid(v).name())

#define warn(cond) {if(cond){ std::cerr << "Warning: (" << # cond << \
                                ") true -- file \"" << __FILE__ << \
                                    "\", line " << __LINE__ << \
                    ", Class '" << DEBUG_TYPE_NAME(*this) <<"'" << std::endl;  }}
 
#define fatal(cond) {if(cond){ std::cerr << "FATAL Error: (" # cond << \
                                ") true -- file\"" << __FILE__ << \
                                  "\", line " << __LINE__ <<  \
                 ", Class '" << DEBUG_TYPE_NAME(*this) << "'" << std::endl; \
                                 abort() ; }}
 
#define WARN(cond) {if(cond){ std::cerr << "Warning: (" << # cond << \
                                ") true -- file \"" << __FILE__ << \
                                    "\", line " << __LINE__ << std::endl; }}
 
#define FATAL(cond) {if(cond){ std::cerr << "FATAL Error: (" # cond << \
                                ") true -- file\"" << __FILE__ << \
                                  "\", line " << __LINE__ << std::endl; \
                                 abort() ; }}


#else  // DEBUG

#define warn(cond)
#define fatal(cond)

#define WARN(cond)
#define FATAL(cond)


#endif // DEBUG

#endif // DEBUG_H
