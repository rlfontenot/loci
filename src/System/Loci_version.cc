#include <Loci_version.h>

namespace Loci {
  
  char *revision_name = "$Name:  $" ;

  std::string version() {
    std::string s = "Loci Library Version: " ;
    std::string rn = revision_name ;
    return s+rn ;
  }
}
