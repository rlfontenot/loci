#ifndef BLACKBOX_H
#define BLACKBOX_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <mpi.h>

#include <Config/conf.h>
#include <Tools/debug.h>
#include <store_rep.h>
#include <istream>
#include <ostream>
#include <data_traits.h>

namespace Loci {
  
  template<class T> class blackboxRepI : public storeRep {
    entitySet store_domain;
    T attrib_data;

  public:
    blackboxRepI() { store_domain = interval(UNIVERSE_MIN,UNIVERSE_MAX); }
    blackboxRepI(const entitySet &p) { store_domain = p;}
    virtual void allocate(const entitySet &p) ;
    virtual ~blackboxRepI();
    virtual store_type RepType() const;
    virtual entitySet domain() const;
    virtual storeRep *new_store(const entitySet &p) const;
    virtual storeRepP remap(const dMap &m) const;
    virtual void copy(storeRepP &st, const entitySet &context);
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context);
    virtual int pack_size(const entitySet &e);
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e);
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual std::ostream &Print(std::ostream &s) const;
    virtual std::istream &Input(std::istream &s);
    virtual void readhdf5(hid_t group, entitySet &en);
    virtual void writehdf5(hid_t group,entitySet& en) const;
    T *get_blackbox() { return &attrib_data; }
  };

  //**************************************************************************/

  template<class T> void blackboxRepI<T>::allocate(const entitySet &p) {
    store_domain = p;
    dispatch_notify();
  }

  //**************************************************************************/

  template<class T> blackboxRepI<T>::~blackboxRepI<T>() {}

  //**************************************************************************/

  template<class T>
    storeRep *blackboxRepI<T>::new_store(const entitySet &p) const
    {
      return new blackboxRepI<T>(p);
    }

  //**************************************************************************/

  template<class T> 
    store_type blackboxRepI<T>::RepType() const 
    {
      return PARAMETER;
    }

  //**************************************************************************/

  template<class T> entitySet blackboxRepI<T>::domain() const {
    return store_domain;
  }

  //**************************************************************************/
        
  template<class T> 
    std::ostream &blackboxRepI<T>::Print(std::ostream &s) const 
    {
      s << '{' << domain() << std::endl;
      Loci::streamoutput(&attrib_data,1,s);
      s << '}' << std::endl;
      return s;
    }

  //**************************************************************************/

  template<class T> 
    std::istream &blackboxRepI<T>::Input(std::istream &s) 
    {
      entitySet e;
      char ch;
    
      do ch = s.get(); while(ch==' ' || ch=='\n');
      if(ch != '{') {
	s.putback(ch);
	e = ~EMPTY;
	allocate(e);
	attrib_data = T();
	Loci::streaminput(&attrib_data,1,s);
	return s;
      }
        
      s >> e;
      allocate(e);
        
      attrib_data = T();
      Loci::streaminput(&attrib_data,1,s);
        
      do ch = s.get(); while(ch==' ' || ch=='\n');
      if(ch != '}') {
	std::cerr << "Incorrect Format while reading blackbox" << std::endl;
	s.putback(ch);
      }

      return s;
    }

  //**************************************************************************/
  template<class T> 
    void blackboxRepI<T>::readhdf5(hid_t group_id, entitySet &user_eset)
    {

    }

  //**************************************************************************/

  template<class T> 
    void blackboxRepI<T>::writehdf5(hid_t group_id, entitySet &eset) const
    {

    }

  //**************************************************************************/

  template<class T> class blackbox : public store_instance {
    typedef blackboxRepI<T> blackboxType;
    T * data;
  public:
    typedef T containerType;
    blackbox() { setRep(new blackboxType); }
    blackbox(blackbox &var) { setRep(var.Rep()); }
    blackbox(storeRepP rp) { setRep(rp); }

    virtual ~blackbox();

    blackbox & operator=(blackbox &p) {setRep(p.Rep()); return *this; }

    blackbox & operator=(storeRepP p) {setRep(p); return *this; }
    blackbox & operator=(const T &v) { *data = v; return *this; }

    virtual void notification();
    
    T * operator->() { return data; }
    const T * operator->() const { return data; }
    
    T * operator&() { return data; }
    const T * operator &() const { return data; }

    T &operator*() { return *data; }
    const T &operator*() const { return *data; }

    T &operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(data == NULL);
      fatal(!Rep()->domain().inSet(indx));
#endif
      return *data;
    }

    const T &operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(data == NULL);
      fatal(!Rep()->domain().inSet(indx));
#endif
      return *data;
    }

    void set_entitySet(const entitySet &ptn) {Rep()->allocate(ptn); }

    entitySet domain() const { return Rep()->domain(); }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s); }
  };

  //**************************************************************************/

  template<class T> blackbox<T>::~blackbox() {}
    
  //**************************************************************************/

  template<class T> 
    void blackbox<T>::notification()
    {  
      NPTR<blackboxType> p(Rep());
      if(p!=0) data = p->get_blackbox();
      warn(p==0);
    }

  //**************************************************************************/

  template<class T> 
    inline std::ostream & operator<<(std::ostream &s, const blackbox<T> &t)
    {
      return t.Print(s);
    }

  //**************************************************************************/

  template<class T> 
    inline std::istream & operator>>(std::istream &s, blackbox<T> &t)
    {
      return t.Input(s);
    }

  //**************************************************************************/

  template<class T> 
    class const_blackbox : public store_instance {
    typedef T containerType;
    typedef blackboxRepI<T> blackboxType;
    const T * data;
    public:
    const_blackbox() { setRep(new blackboxType); }
    const_blackbox(const_blackbox<T> &var) { setRep(var.Rep()); }
    const_blackbox(blackbox<T> &var) { setRep(var.Rep()); }
    const_blackbox(storeRepP rp) { setRep(rp); }
    
    virtual ~const_blackbox();

    const_blackbox & operator=(const_blackbox<T> &p)
    { setRep(p.Rep); return *this;}
    const_blackbox & operator=(blackbox<T> &p)
    { setRep(p.Rep); return *this;}
    const_blackbox & operator=(storeRepP p)
    { setRep(p); return *this;}

    virtual void notification();
    virtual instance_type access() const;
        
    const T * operator->() const { return data; }
    
    const T * operator &() const { return data; }

    const T &operator*() const { return *data; }

    const T &operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(data == NULL);
      fatal(!Rep()->domain().inSet(indx));
#endif
      return *data;
    }

    entitySet domain() const { return Rep()->domain(); }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
  };

  //**************************************************************************/

  template<class T> const_blackbox<T>::~const_blackbox() {}

  //**************************************************************************/

  template<class T> 
    void const_blackbox<T>::notification() 
    {  
      NPTR<blackboxType> p(Rep());
      if(p!=0) data = p->get_blackbox();
      warn(p==0);
    }
    
  //**************************************************************************/

  template<class T> 
    storeRepP blackboxRepI<T>::remap(const dMap &m) const 
    {
      blackbox<T> r;
      r.set_entitySet(m.image(m.domain()&domain()));
      *r = attrib_data;
      return r.Rep();
    }

  //**************************************************************************/

  template<class T> 
    void blackboxRepI<T>::copy(storeRepP &st, const entitySet &context) 
    {
      blackbox<T> p(st);
      attrib_data = *p;
      warn((store_domain - context) != EMPTY);
      store_domain = context;
      dispatch_notify();
    }

  //**************************************************************************/

  template<class T> 
    void blackboxRepI<T>::gather(const dMap &m, storeRepP &st,
				 const entitySet &context) 
    {

    }

  //**************************************************************************/

  template<class T> 
    void blackboxRepI<T>::scatter(const dMap &m, storeRepP &st,
				  const entitySet &context) 
    {

    }

  //**************************************************************************/
 
  template <class T> 
    int blackboxRepI<T>::pack_size(const entitySet &eset) 
    {
      return 0;
    }

  //**************************************************************************/

/*   template <class T> */
/*     int blackboxRepI<T>::get_mpi_size(IDENTITY_CONVERTER c, */
/* 				      const entitySet &eset) */
/*     { */

/*     } */

  //**************************************************************************/

/*   template <class T> */
/*     int blackboxRepI<T>::get_mpi_size(USER_DEFINED_CONVERTER c, */
/* 				      const entitySet &eset) */
/*     { */

/*     } */
  //**************************************************************************/

  template <class T> 
    void blackboxRepI<T>::pack(void *ptr, int &loc, int &size,
			       const entitySet &e) 
    {

    }

  //**************************************************************************/
/*   template <class T> */
/*     void blackboxRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf, */
/* 				    int &position, int outcount) */
/*     { */

/*     } */

  //**************************************************************************/

/*   template <class T>  */
/*     void blackboxRepI<T>::packdata(USER_DEFINED_CONVERTER c, void *outbuf,  */
/* 				    int &position, int outcount) */
/*     { */

/*     } */
  //**************************************************************************/

  template <class T> 
    void blackboxRepI<T>::unpack(void *ptr, int &loc, int &size,
				 const sequence &seq)
    {

    }  


  //**************************************************************************/
/*   template <class T>  */
/*     void blackboxRepI<T>::unpackdata(IDENTITY_CONVERTER c, void *inbuf, */
/* 				      int &position, int &insize) */
/*     { */

/*     } */

  //***********************************************************************/
/*   template <class T>  */
/*     void blackboxRepI<T>::unpackdata(USER_DEFINED_CONVERTER c, void *inbuf,  */
/* 				      int &position, int &insize) */
/*     { */

/*     } */

  //***********************************************************************/
  
  template<class T>
    store_instance::instance_type const_blackbox<T>::access() const
    {
      return READ_ONLY;
    }

  //**************************************************************************/
    
  template<class T> 
    inline std::ostream & operator<<(std::ostream &s,
				     const const_blackbox<T> &t)
    {
      return t.Print(s);
    }

  //**************************************************************************/
/*   template <class T>  */
/*     void blackboxRepI<T> :: hdf5write(hid_t group_id, IDENTITY_CONVERTER g, */
/* 				       const entitySet &eset) const */
/*     { */

/*     } */

  //*********************************************************************/

/*   template <class T>  */
/*     void blackboxRepI<T> :: hdf5write(hid_t group_id, USER_DEFINED_CONVERTER g,  */
/* 				       const entitySet &eset) const */
/*     { */

/*     } */

  //**************************************************************************/
/*   template <class T>  */
/*     void blackboxRepI<T> :: hdf5read(hid_t group_id, IDENTITY_CONVERTER g) */
/*     {  */

/*     } */

  //*************************************************************************/

/*   template <class T>  */
/*     void blackboxRepI<T> :: hdf5read(hid_t group_id, USER_DEFINED_CONVERTER g) */
/*     {  */

/*     } */

  //***************************************************************************

}

#endif
    
