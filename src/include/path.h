#ifndef PATH_H
#define PATH_H

#include <Tools/debug.h>
#include <store_rep.h>

namespace Loci {
    
    class pathRep : public storeRep {
        entitySet store_domain ;
        sequence path ;
      public:
        pathRep() ;
        pathRep(const sequence &p) ;
        virtual ~pathRep() ;
        virtual void allocate(const entitySet &p) ;
        virtual storeRep *new_store(const entitySet &p) const ;
        virtual store_type RepType() const ;
        virtual const entitySet &domain() const ;
        virtual std::ostream &Print(std::ostream &s) const ;
        virtual std::istream &Input(std::istream &s) ;
        sequence *get_path() { return &path ; }
    } ;

    class path : public store_instance {
        typedef pathRep pathType ;
        sequence *data ;
      public:
        path() ;
        path(path &var) ;
        path(const sequence &ptn) ;
        path(storeRepP rep) { setRep(rep) ;}
        virtual ~path() ;

        path & operator=(path &p)
        { setRep(p.Rep()) ; return *this ;}
        path & operator=(storeRepP p)
        { setRep(p) ;  return *this ;}
        path & operator=(const sequence &v)
        { *data = v ; return *this ; }

        virtual void notification() ;
    
        sequence * operator&() { return data ; }
        const sequence * operator &() const { return data ; }

        sequence *operator->() { return data ; }
        const sequence *operator->() const { return data ; }
        
        sequence &operator*() { return *data ; }
    
        operator storeRepP() { return Rep() ; }

        operator sequence() { return *data ; }
        operator sequence() const { return *data ; }
        std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
        std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
    } ;

    inline std::ostream & operator<<(std::ostream &s, const path &t)
        { return t.Print(s) ; }

    inline std::istream & operator>>(std::istream &s, path &t)
        { return t.Input(s) ; }

}

#endif
    
