#include <constraint.h>
#include <Tools/stream.h>
#include <Map.h>
#include <DMap.h>

namespace Loci {



  constraintRep::constraintRep()
  {
  }

  constraintRep::constraintRep(const entitySet &p)
  {
    constraint_set = p ;
  }

  constraintRep::~constraintRep()
  {
  }

  void constraintRep::allocate(const entitySet &p)
  {
    constraint_set = p ;
    dispatch_notify() ;
  }

  storeRep *constraintRep::new_store(const entitySet &p) const
  {
    return new constraintRep(p) ;
  }
  storeRep *constraintRep::new_store(const entitySet &p, const int* cnt) const
  {
    storeRep* sp ;
    cerr << " This method should not be called for a constraint " << endl ;
    return sp ;
  }

  storeRepP constraintRep::remap(const dMap &m) const {
    entitySet newconstraint = m.image(m.domain()&constraint_set) ;
    constraint r ;
    r = newconstraint ;
    return r.Rep() ;
  }

  void constraintRep::copy(storeRepP &st, const entitySet &context) {
    constraint cs(st) ;
    entitySet sent,tent ;
    sent = *cs ;
    tent = constraint_set ;
    tent -= context ;
    tent += sent & context ;
    constraint_set = tent ;
    dispatch_notify() ;
  }

  void constraintRep::gather(const dMap &m, storeRepP &st,
                             const entitySet &context) {
    constraint cs(st) ;
    entitySet tent = constraint_set ;
    tent -= context ;
    entitySet img = *cs ;
    FORALL(context,i) {
      if(img.inSet(m[i]))
        tent+=i ;
    } ENDFORALL ;
    constraint_set = tent ;
    dispatch_notify() ;
  }

  void constraintRep::scatter(const dMap &m, storeRepP &st,
                              const entitySet &context) {
    constraint cs(st) ;
    entitySet map_image = m.image(context) ;
    entitySet tent = constraint_set ;
    tent -= map_image ;
    entitySet img = *cs ;
    tent += m.image(context&img) ;
    constraint_set = tent ;
    dispatch_notify() ;
  }
  
  int constraintRep::pack_size(const entitySet &e) {
    warn(true) ;
    return 0 ;
  }
  
  void constraintRep::pack(void *ptr, int &loc, int &size, const entitySet&e) {
    warn(true) ;
  }
  
  void constraintRep::unpack(void *ptr, int &loc, int &size, const sequence &seq) {
    warn(true);
  }
  
  store_type constraintRep::RepType() const
  {
    return CONSTRAINT ;
  }

  entitySet constraintRep::domain() const {
    return constraint_set ;
  }

  ostream &constraintRep::Print(ostream &s) const {
    s << constraint_set << endl ;
    return s ;
  }
  DatatypeP constraintRep::getType() {
    return DatatypeP(new AtomicType(INT)) ;
  }
  frame_info constraintRep::read_frame_info(hid_t group_id) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  frame_info constraintRep::write_frame_info(hid_t group_id) {
    warn(true) ; 
    frame_info fi ;
    return fi ;
  }
  
  void constraintRep::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en){
    warn(true) ;
    
/*
    try{
      //get constraint data
      H5::DataSet dataset_constraint = group.openDataSet( "constraint");
      H5::DataSpace dataspace_constraint = dataset_constraint.getSpace();
      hsize_t dims_constraint[1];
      dataspace_constraint.getSimpleExtentDims( dims_constraint, NULL);
      int *data_constraint = new int[dims_constraint[0]];
      dataset_constraint.read( data_constraint, H5::PredType::NATIVE_INT );
      for(int i=0;i<dims_constraint[0];i++){
        constraint_set |=interval(data_constraint[i],data_constraint[i+1]);
        i++;
      }
      delete [] data_constraint;
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
*/

  }

  void constraintRep::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const{
    warn(true) ;
    /*
    hsize_t dimf_constraint[1];
    int RANK=1;
    
    int num_intervals=constraint_set.num_intervals();
    dimf_constraint[0]=num_intervals*2;
    interval *it = new interval[num_intervals];
    int *data_constraint = new int[num_intervals*2];//get the constraint data
    for(int i=0;i<num_intervals;i++){
      it[i]=constraint_set[i];
      data_constraint[i*2]=it[i].first;
      data_constraint[i*2+1]=it[i].second;
    }
    try{
      //write constraint
      H5::DataSpace dataspace_constraint( RANK, dimf_constraint );
      H5::DataSet dataset_constraint = group.createDataSet( "constraint", H5::PredType::NATIVE_INT, dataspace_constraint );
      dataset_constraint.write( data_constraint, H5::PredType::NATIVE_INT );
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
    
    delete [] it;
    delete [] data_constraint;
*/
  }

  istream &constraintRep::Input(istream &s) {
    entitySet e ;
    s >> e ;
    allocate(e) ;
    return s ;
  }

  constraint::constraint()
  {
    setRep(new constraintType) ;
  }
    
  constraint::constraint(const constraint &var)
  {
    setRep(var.Rep()) ;
  }

  constraint::~constraint()
  {
  }

  void constraint::notification()
  {
    NPTR<constraintType> p(Rep());
    if(p!=0)
      data = p->get_constraint() ;
    warn(p==0);
  }

}
