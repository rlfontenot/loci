#ifndef HDF5_MEMENTO_H_
#define HDF5_MEMENTO_H_

#include <typeinfo>
#include <hdf5_traits.h>
#include <hdf5CC/H5cpp.h>

#include <vector>
#include <algorithm>


namespace Loci {

template<class T>
class Memento
{ };


template< class T>
class Memento< std::vector<T> >
{
   typedef hdf5_schema_converter_traits<std::vector<T> > converter_traits; 

public:
  explicit Memento( const std::vector<T> &buf)
  { obj = buf; };

  int  getSize()  const
  { return obj.size(); };

  void getState(typename converter_traits::memento_type *buf, int &size)
  {
    size = obj.size();
    for( int i = 0; i < size; i++) 
         buf[i] = obj[i];
  }
  std::vector<T> setState(typename converter_traits::memento_type *buf, int size)
  {
    obj.clear();
 
    for( int i = 0; i < size; i++)
         obj.push_back( buf[i] );

    return( obj );
  }

private:
   std::vector<T>  obj;
};

















//*************************************************************************/
// set store memento
//*************************************************************************/

/*
template <class T>
inline void set_store_memento( std::vector<T> *attrib_data, 
                               const entitySet en,
                               Memento< std::vector<T> >& memento)
{
  cout << " Commented for the time being " << endl;
  Memento< std::vector<T> >::var_type* temp_data = memento.get_variable();
  int* temp_framing = memento.get_framing();

  //restore data in original object
  int gi=0,ti=0;
  
  FORALL(en,ii){
    for(int g=0;g<temp_framing[gi];g++){
      attrib_data[ii].push_back(temp_data[ti]);
      ti++;
    }
    gi++;
  }ENDFORALL
};
*/

//****************************************************************************/

/*
template <class T>
inline void set_store_memento( std::complex<T>* attrib_data, 
                               const entitySet en,
                               Memento< std::complex<T> >& memento)
{
  cout << " Commented for the time being " << endl;

  Memento<std::complex<T> >::var_type* temp_data = memento.get_variable();
  int* temp_framing = memento.get_framing();

  //restore data in original object
  int gi=0,ti=0;
  FORALL(en,ii){
      for(int g=0;g<temp_framing[gi];g++){
        std::complex<T> x(temp_data[ti],temp_data[ti++]);
        attrib_data[ii] = x;
        ti++;
      }
      gi++;
 } ENDFORALL
};

*/

//***************************************************************************/
//create storeVec memento
//***************************************************************************/

/*
template <class W>
void create_storeVec_memento(const W* attrib_data,const entitySet en,int size,
                                                Memento<W>& memento)
{
    cerr<<"create_storeVec_memento<class T> should not be called"<<endl;
}
*/

//****************************************************************************/

/*
template <class T>
inline void create_storeVec_memento( const std::vector<T> *attrib_data,
                                     const entitySet en, int vecsize, 
									          Memento< std::vector<T> >& memento)
{
  std::vector<int>  veclength;
  std::vector<T>    vec;

  entitySet:: const_iterator  ci;

  int indx;
  int arraySize = 0, numBuckets = 0;
  for( ci = en.begin(); ci != en.end(); ++ci) {
       numBuckets += vecsize;
       for( int i = 0; i < vecsize; i++) {
            indx       = (*ci)*vecsize + i;
            vec        = attrib_data[indx];
            arraySize += vec.size();
            veclength.push_back( vec.size() );
       }
  }

//----------------------------------------------------------------------------
// Create an array of which contains the size of each vector
//----------------------------------------------------------------------------
 
  int  *veclen = new int[numBuckets];
  for( int i = 0; i < numBuckets; i++)
       veclen[i] = veclength[i];

  veclength.clear();


//----------------------------------------------------------------------------
// Store the data in 1D
//----------------------------------------------------------------------------
  Memento< std::vector<T> >::var_type *data;
 
  data = new Memento< std::vector<T> >::var_type[arraySize];

  int k = 0;
  for( ci = en.begin(); ci != en.end(); ++ci) {
       for( int i = 0; i < vecsize; i++) {
            indx  = (*ci)*vecsize + i;
            vec   = attrib_data[indx];
            for( int j = 0; j < vec.size(); j++)
                 data[k++] = vec[j];
       }
  }

  memento.setState(data, arraySize);
  memento.setBucket( veclen, numBuckets );

};
*/

//****************************************************************************/
 
/*
template <class T>
inline void create_storeVec_memento(const std::complex<T>* attrib_data, 
                                    const entitySet en, int vecsize,
                                    Memento< std::complex<T> >& memento)
{

  int arraySize = 2*en.size()*vecsize;

  Memento< std::complex<T> >::var_type *data;
 
  data = new Memento< std::complex<T> >::var_type[arraySize];

  int k = 0, indx;
  entitySet:: const_iterator  ci;
  for( ci = en.begin(); ci != en.end(); ++ci) {
       for( int i = 0; i < vecsize; i++) {
            indx  = (*ci)*vecsize + i;
            data[k++] = attrib_data[indx].real();
            data[k++] = attrib_data[indx].imag();
       }
  }
  memento.setState(data, arraySize);

  int numBuckets = 1;
  int veclen     = 2;
  
  memento.setBucket( &veclen, numBuckets );
};
*/

//***************************************************************************/
// set storeVec memento
//***************************************************************************/

/*
template <class W>
void set_storeVec_memento( W* attrib_data, const entitySet en, int size, 
                           Memento<W>& memento)
{
    cerr<<"set_storeVec_memento<class T> should not be called"<<endl;
};
*/

//***************************************************************************/

/*
template <class T>
inline void set_storeVec_memento( std::vector<T>* attrib_data, 
                                  const entitySet en, int size,
                                  Memento< std::vector<T> >& memento)
{
    Memento< std::vector<T> >::var_type* temp_data = memento.get_variable();
   int* temp_framing = memento.get_framing();

   //restore data in original object
   int gi=0,ti=0;
   FORALL(en,ii){
     for(int j=0;j<size;j++){
       for(int g=0;g<temp_framing[gi];g++){
	 attrib_data[ii*size+j].push_back(temp_data[ti]);
	 ti++;
       }
       gi++;
     }
   }ENDFORALL
};

*/

//***************************************************************************/

/*
template <class T>
inline void set_storeVec_memento( std::complex<T>* attrib_data, 
                                  const entitySet en, int size,
                                  Memento< std::complex<T> >& memento)
{
    Memento< std::complex<T> >::var_type* temp_data = memento.get_variable();
   int* temp_framing = memento.get_framing();

   //restore data in original object
   int gi=0,ti=0;
   FORALL(en,ii){
     for(int j=0;j<size;j++){
       for(int g=0;g<temp_framing[gi];g++){
	 std::complex<T> x(temp_data[ti],temp_data[ti++]);
	 attrib_data[ii*size+j]= x;
	 ti++;
       }
       gi++;
     }
   }ENDFORALL
};

*/

//***************************************************************************/
// set param memento
//***************************************************************************/

/*
template <class W>
void set_param_memento( W& attrib_data, const entitySet en, 
                        Memento<W>& memento)
{
   cerr<<"set_param_memento<class T> should not be called"<<endl;
};
*/

//***************************************************************************/

/*
template <class T>
inline void set_param_memento( std::vector<T>& attrib_data, 
                               const entitySet en,
                               Memento< std::vector<T> >& memento)
{
  Memento< std::vector<T> >::var_type* temp_data = memento.get_variable();
   int* temp_framing = memento.get_framing();

   //restore data in original object
   int gi=0;
   for(int g=0;g<temp_framing[gi];g++){
     attrib_data.push_back(temp_data[g]);
   }
};
*/

//***************************************************************************/
 
/*
template <class T>
inline void set_param_memento( std::complex<T>& attrib_data,
                               const entitySet en,
                               Memento< std::complex<T> >& memento)
{
   Memento< std::complex<T> >::var_type* temp_data = memento.get_variable();
   int* temp_framing = memento.get_framing();

   //restore data in original object
   int gi=0;
   for(int g=0;g<temp_framing[gi];g++){
     std::complex<T> x(temp_data[gi],temp_data[gi++]);
     attrib_data = x;
     gi++;
   }
};
*/

//***************************************************************************/
// memento hdf5write
//***************************************************************************/

template <class T, class W>
void memento_hdf5write(H5::Group group, T t, const Memento<W> &memento, 
                       const entitySet &en, int size)
{
   cerr<<" memento_hdf5write:This should not be called!"<<endl;
};

//***************************************************************************/

template <class W>
void memento_hdf5write(H5::Group group, Need_Variable_data t, 
                       const Memento<W> &memento, const entitySet &en, 
					        int size)
{
   cout << " Calling Useless function " << endl;
   exit(0);
}

//***************************************************************************/

template <class W>
void memento_hdf5write(H5::Group group, Memento<W> &memento )
{

  typedef    hdf5_schema_converter_traits <W> _schema_traits_type;

  hsize_t    dimension[1];

  int rank   = 1; 

  dimension[0] = memento.getSize();

  try {
      //-----------------------------------------------------------------------
      // Write (variable) Data ....
      //-----------------------------------------------------------------------

      H5::DataSpace vDataspace( rank, dimension );
      H5::DataType  vDatatype = _schema_traits_type::get_variable_HDF5_type();
      H5::DataSet   vDataset  = group.createDataSet( "variable", vDatatype, vDataspace);

      vDataset.write( memento.getState(), vDatatype );

      //----------------------------------------------------------------------
      // Write plane data type information
      //----------------------------------------------------------------------
      /*
      dimension[0] = 1;
      H5::DataSpace POD_dataspace( rank, dimension );
      H5::DataType  POD_datatype = _schema_traits_type::get_POD_HDF5_type();
      H5::DataSet   POD_dataset  = group.createDataSet( "POD", POD_datatype, POD_dataspace);

      Memento<W>::info_type pod_info = memento.get_info();
      POD_dataset.write( &temp_info, POD_datatype );
      */

      //-----------------------------------------------------------------------
      // Write Component Size 
      //-----------------------------------------------------------------------
/*

      dimension[0]   = memento.numBuckets();

      H5::DataType  frameDatatype = H5::PredType::NATIVE_INT;
      H5::DataSpace frameDataspace( rank, dimension );
      H5::DataSet   frameDataset = group.createDataSet( "framing", 
                                                        frameDatatype, frameDataspace);
      frameDataset.write( memento.getBucket(), frameDatatype );
*/

      //-----------------------------------------------------------------------
      // HDF5 write over 
      //-----------------------------------------------------------------------

   }
   catch( H5::HDF5DatasetInterfaceException error   ) { error.printerror(); }
   catch( H5::HDF5DataspaceInterfaceException error ) { error.printerror(); }
   catch( H5::HDF5DatatypeInterfaceException error  ) { error.printerror(); }

 } 

//***************************************************************************/

template <class W> 
void memento_hdf5write(H5::Group group, No_Need_Variable_data t, 
                       const Memento<W> &memento, 
                       const entitySet &en, int size)
{
/*
   typedef hdf5_schema_converter_traits <W> _schema_traits_type;
   hsize_t dimf_POD[1];

   int rank     = 1;
   dimf_POD[0]  = 1;

   Memento<W>::info_type temp_info = memento.get_info();

   //write data,only POD data is needed to write out    
   try
   {
      H5::DataType POD_datatype=_schema_traits_type::get_POD_HDF5_type();
      H5::DataSpace POD_dataspace( rank, dimf_POD );
      H5::DataSet POD_dataset = group.createDataSet( "POD", POD_datatype, POD_dataspace);
      POD_dataset.write( &temp_info, POD_datatype );
   }

   catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
   catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
   catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
*/
};

//***************************************************************************/
// memento hdf5read
//***************************************************************************/

template <class W>
void memento_hdf5read(H5::Group group, Need_Variable_data t, 
                      Memento<W> &memento, entitySet &en)
{
/*
  typedef hdf5_schema_converter_traits <W> _schema_traits_type;
  
    try{
	//get POD data
        H5::DataType datatype_POD=_schema_traits_type::get_POD_HDF5_type();
	H5::DataSet dataset_POD = group.openDataSet( "POD");
	H5::DataSpace dataspace_POD = dataset_POD.getSpace();
	hsize_t dims_POD[1];
	dataspace_POD.getSimpleExtentDims( dims_POD, NULL);
	Memento<W>::info_type *data_POD = new Memento<W>::info_type[dims_POD[0]];
	dataset_POD.read( data_POD, datatype_POD );
	memento.set_info(*data_POD);

	//get variable data
        H5::DataType datatype_variable=_schema_traits_type::get_variable_HDF5_type();
	H5::DataSet dataset_variable = group.openDataSet( "variable");
	H5::DataSpace dataspace_variable = dataset_variable.getSpace();
	hsize_t dims_variable[1];
	dataspace_variable.getSimpleExtentDims( dims_variable, NULL);
	Memento<W>::var_type *data_variable = new Memento<W>::var_type[dims_variable[0]];
	dataset_variable.read( data_variable, datatype_variable);
	memento.set_variable(data_variable);

	//get framing data
        H5::DataType datatype_framing = H5::PredType::NATIVE_INT;
	H5::DataSet dataset_framing = group.openDataSet( "framing");
	H5::DataSpace dataspace_framing = dataset_framing.getSpace();
	hsize_t dims_framing[1];
	dataspace_framing.getSimpleExtentDims( dims_framing, NULL);
	int *data_framing = new int[dims_framing[0]];
	dataset_framing.read( data_framing, H5::PredType::NATIVE_INT );
	memento.set_framing(data_framing);

	delete [] data_POD;
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
*/
};

//***************************************************************************/

template <class W>
void memento_hdf5read( H5::Group group, No_Need_Variable_data t, 
                       Memento<W> &memento, entitySet &en)
{
/*
    typedef hdf5_schema_converter_traits <W> _schema_traits_type;
    try{
	//get POD data
        H5::DataType datatype_POD=_schema_traits_type::get_POD_HDF5_type();
	H5::DataSet dataset_POD = group.openDataSet( "POD");
	H5::DataSpace dataspace_POD = dataset_POD.getSpace();
	hsize_t dims_POD[1];
	dataspace_POD.getSimpleExtentDims( dims_POD, NULL);
	Memento<W>::info_type *data_POD = new Memento<W>::info_type[dims_POD[0]];
	dataset_POD.read( data_POD, datatype_POD );
	memento.set_info(*data_POD);

	delete [] data_POD;;
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
*/
};

//***************************************************************************/

}
#endif
