#include <constraint.h>
#include <Tools/stream.h>

namespace Loci {



constraintRep::constraintRep()
{
}

constraintRep::constraintRep(const entitySet &p)
{
    constraint = p ;
}

constraintRep::~constraintRep()
{
}

void constraintRep::allocate(const entitySet &p)
{
    constraint = p ;
    dispatch_notify() ;
}

storeRep *constraintRep::new_store(const entitySet &p) const
{
    return new constraintRep(p) ;
}

store_type constraintRep::RepType() const
{
    return CONSTRAINT ;
}

const entitySet &constraintRep::domain() const {
    return constraint ;
}

ostream &constraintRep::Print(ostream &s) const {
    s << constraint << endl ;
    return s ;
}

void constraintRep::readhdf5( H5::Group group){
  try{
     //get constraint data
     H5::DataSet dataset_constraint = group.openDataSet( "constraint");
     H5::DataSpace dataspace_constraint = dataset_constraint.getSpace();
     hsize_t dims_constraint[1];
     dataspace_constraint.getSimpleExtentDims( dims_constraint, NULL);
     int *data_constraint = new int[dims_constraint[0]];
     dataset_constraint.read( data_constraint, H5::PredType::NATIVE_INT );
     for(int i=0;i<dims_constraint[0];i++){
       constraint |=interval(data_constraint[i],data_constraint[i+1]);
       i++;
     }
     delete [] data_constraint;
  }
  catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
  catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
  catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
}

void constraintRep::writehdf5( H5::Group group,entitySet& en) const{
  hsize_t dimf_constraint[1];
   int RANK=1;
    
    int num_intervals=constraint.num_intervals();
    dimf_constraint[0]=num_intervals*2;
    interval *it = new interval[num_intervals];
    int *data_constraint = new int[num_intervals*2];//get the constraint data
    for(int i=0;i<num_intervals;i++){
      it[i]=constraint[i];
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
    
constraint::constraint(constraint &var)
{
    setRep(var.Rep()) ;
}

constraint::constraint(const entitySet &ptn)
{
    setRep(new constraintType(ptn)) ;
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
