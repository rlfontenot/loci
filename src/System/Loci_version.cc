#include <Loci_version.h>

namespace Loci {
  
  char *revision_name = "$Name: test-1-1e $" ;

  std::string version() {
    char *p = revision_name;
    while(*p!=':' && *p!='\0')
      ++p ;
    if(*p!= '\0')
      ++p ;
    while(*p!=' ' && *p!='\0')
      ++p ;
    if(*p!= '\0')
      ++p ;
    std::string rn ;
    while(*p!='$' &&  *p!=' ' && *p!='\0') 
      rn += *p++ ;
    return rn ;
  }
}
