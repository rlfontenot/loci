#include <Loci_version.h>

namespace Loci {
  
  const char *revision_name = "$Name:  $" ;

  std::string version() {
    const char *p = revision_name;
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
