#ifndef STORE_H
#define STORE_H 1

#include <Tools/debug.h>
#include <store_rep.h>

namespace Loci {

    template<class T> class storeRepI : public storeRep {
        T *alloc_pointer ;
        T *base_ptr ;
        entitySet store_domain ;
      public:
        storeRepI() { alloc_pointer = 0 ; base_ptr = 0 ; }
        storeRepI(const entitySet &p) { alloc_pointer=0 ; allocate(p) ; }
        virtual void allocate(const entitySet &ptn) ;
        virtual ~storeRepI()  ;
        virtual storeRep *new_store(const entitySet &p) const ;
        virtual store_type RepType() const ;
        virtual std::ostream &Print(std::ostream &s) const ;
        virtual std::istream &Input(std::istream &s) ;
        virtual const entitySet &domain() const ;
        T * get_base_ptr() const { return base_ptr ; }
    } ;

    template<class T> void storeRepI<T>::allocate(const entitySet &ptn) {
        if(alloc_pointer) delete[] alloc_pointer ;
        alloc_pointer = 0 ;
        base_ptr = 0 ;
        if(ptn != EMPTY) {
            int top = ptn.Min() ; int size = ptn.Max()-top+1 ;
            alloc_pointer = new(T[size]) ;
            base_ptr = alloc_pointer - top ;
        }
        store_domain = ptn ;
        dispatch_notify() ;
    }

    template<class T> std::ostream &storeRepI<T>::Print(std::ostream &s) const {
        s << '{' << domain() << std::endl ;
        FORALL(domain(),ii) {
            s << base_ptr[ii] << std::endl ;
        }ENDFORALL ;
        s << '}' << std::endl ;
        return s ;
    }
    


    template<class T> std::istream &storeRepI<T>::Input(std::istream &s) {
        entitySet e ;
        char ch ;
        
        do ch = s.get(); while(ch==' ' || ch=='\n') ;
        if(ch != '{') {
            cerr << "Incorrect Format while reading store" << std::endl ;
            s.putback(ch) ;
            return s ;
        }
        s >> e ;
        allocate(e) ;
        
        FORALL(e,ii) {
            s >> base_ptr[ii] ;
        } ENDFORALL ;
        
        do ch = s.get(); while(ch==' ' || ch=='\n') ;
        if(ch != '}') {
            cerr << "Incorrect Format while reading store" << std::endl ;
            s.putback(ch) ;
        }
        return s ;
    }

    template<class T>  storeRepI<T>::~storeRepI<T>() {
        if(alloc_pointer) delete[] alloc_pointer ;
    }
    
    template<class T>  const entitySet &storeRepI<T>::domain() const {
        return store_domain ;
    }

    template<class T>
        storeRep *storeRepI<T>::new_store(const entitySet &p) const {
        return new storeRepI<T>(p)  ;
    }

    template<class T> store_type storeRepI<T>::RepType() const {
        return STORE ;
    }

    template<class T> class store : public store_instance {
        typedef storeRepI<T> storeType ;
        T* base_ptr ;
      public:
        store() { setRep(new storeType); }
        store(store &var) { setRep(var.Rep()) ; }
        store(const entitySet &ptn) {
            setRep(new storeType(ptn)) ;
        }

        virtual ~store() ;
        virtual void notification() ;

        store<T> & operator=(store<T> &str) {
            setRep(str.Rep()) ;
            return *this ;
        }

        store<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

        void initialize(const entitySet &ptn) { Rep()->allocate(ptn) ; }
        void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

        entitySet domain() const { return Rep()->domain(); }
        std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
        std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

        operator storeRepP() { return Rep() ; }
        T &elem(int indx) {
#ifdef BOUNDS_CHECK
            fatal(base_ptr==NULL); 
            fatal(!Rep()->domain().inSet(indx)) ;
#endif 
            return base_ptr[indx]; }
        const T &elem(int indx) const {
#ifdef BOUNDS_CHECK
            fatal(base_ptr==NULL); 
            fatal(!Rep()->domain().inSet(indx)) ;
#endif 
            return base_ptr[indx]; }
            
        T &operator[](int indx) { return elem(indx); }
        const T&operator[](int indx) const { return elem(indx); }

    } ;

    template<class T> store<T>::~store<T>() { }
    
    template<class T> void store<T>::notification() {
        NPTR<storeType> p(Rep()) ;
        if(p != 0)
          base_ptr = p->get_base_ptr() ;
        warn(p == 0) ;
    }

    template<class T> inline std::ostream & operator<<(std::ostream &s, const store<T> &t)
        { return t.Print(s) ; }

    template<class T> inline std::istream & operator>>(std::istream &s, store<T> &t)
        { return t.Input(s) ; }



    template<class T> class const_store : public store_instance {
        typedef storeRepI<T> storeType ;
        const T * base_ptr ;
      public:
        const_store() { setRep(new storeType) ; }
        const_store(store<T> &var) { setRep(var.Rep()) ; }
        const_store(const_store &var) { setRep(var.Rep()) ; }
        const_store(const entitySet &ptn) { setRep(new storeType(ptn)) ; }

        virtual ~const_store() ;
        virtual void notification() ;


        virtual instance_type access() const  ;
        
        const_store<T> & operator=(const_store<T> &str) {
            setRep(str.Rep()) ;
            return *this ;
        }

        const_store<T> & operator=(store<T> &str) {
            setRep(str.Rep()) ;
            return *this ;
        }

        const_store<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

        entitySet domain() const { return Rep()->domain(); }
        std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }

        operator storeRepP() { return Rep() ; }

        const T &elem(int indx) const {
#ifdef BOUNDS_CHECK
            fatal(base_ptr==NULL); 
            fatal(!Rep()->domain().inSet(indx)) ;
#endif 
            return base_ptr[indx]; }
            
        const T&operator[](int indx) const { return elem(indx); }

    } ;

    template<class T> const_store<T>::~const_store<T>() { }
    
    template<class T> void const_store<T>::notification() {
        NPTR<storeType> p(Rep()) ;
        if(p != 0)
          base_ptr = p->get_base_ptr() ;
        warn(p == 0) ;
    }

    template<class T> store_instance::instance_type
        const_store<T>::access() const { return READ_ONLY; }
        

}

#endif
