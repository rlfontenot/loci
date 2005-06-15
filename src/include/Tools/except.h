#ifndef EXCEPT_H
#define EXCEPT_H

#include <string>
#include <iostream>

namespace Loci {

    class BasicException {
      public:
        int code ;
        BasicException() {code = 0 ;}
        BasicException(int val) : code(val) {}
        virtual ~BasicException() {}
        virtual std::ostream &Print(std::ostream &s) const {
            s << "Loci exception, code = " << code << std::endl ;
            return s ;
        }
    } ;

    struct StringError : public BasicException {
      public:
        std::string msg ;
        StringError() {}
        StringError(std::string m) :msg(m) {code = 0 ;}
        virtual ~StringError() {}
        virtual std::ostream &Print(std::ostream &s) const {
            s << msg << std::endl ;
            return s ;
        }
    } ;

    struct StringWarning : public BasicException {
      public:
        std::string msg ;
        StringWarning() {}
        StringWarning(std::string m) :msg(m) {code = 0 ;}
        virtual std::ostream &Print(std::ostream &s) const {
            s << msg << std::endl ;
            return s ;
        }
    } ;

}

#endif

    
