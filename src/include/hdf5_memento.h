#ifndef HDF5_MEMENTO_H_
#define HDF5_MEMENTO_H_

#include <hdf5_traits.h>
#include <hdf5CC/H5cpp.h>
#include <Tools/stream.h>
#include <typeinfo>
namespace Loci {

 template <class W> class Memento {
  public:
    typedef hdf5_schema_converter_traits<W> schema_converter;
    typedef typename schema_converter::POD_memento_type   info_type;
    typedef typename schema_converter::variable_memento_type   variable_type;
    void set_info(info_type info){info_data=info;};
    void set_variable(variable_type* var_data){variable_data=var_data;};
    void set_framing(int* fm){framing=fm;};
    info_type get_info() const {return info_data;};
    variable_type* get_variable() const {return variable_data;};
    int* get_framing() const {return framing;};
  private:
    info_type info_data;
    variable_type* variable_data;
    int* framing;
  };

  //-------------------create store memento-------------------------//
  //----------------------------------------------------------//
  template <class W> void create_store_memento(const W* original_obj,const entitySet en, Memento<W>& memento){
    cerr<<"create_store_memento<class T> should not be called"<<endl;
  }

  template <class T> inline void create_store_memento(const std::vector<T>* original_obj,const entitySet en, Memento< std::vector<T> >& memento){
      std::vector<T> temp_data;
      int gi=0;
      int* temp_framing =new int[en.size()] ;
      Memento< std::vector<T> >::info_type temp_info=0;
      FORALL(en,ii){
	temp_data=original_obj[ii];
	temp_info += temp_data.size();
	temp_framing[gi]=temp_data.size();
	gi++;
      }ENDFORALL
	 memento.set_info(temp_info);
	 memento.set_framing(temp_framing);
	 Memento< std::vector<T> >::variable_type*  temp_variable= new Memento< std::vector<T> >::variable_type[memento.get_info()];
      int ti=0;
      gi=0;
      FORALL(en,ii){
	temp_data=original_obj[ii];
	for(int g=0;g<temp_framing[gi];g++){
	  temp_variable[ti]=temp_data[g];
	  ti++;
	}
	gi++;
      }ENDFORALL
	 memento.set_variable(temp_variable);

  };
  
  template <class T> inline void create_store_memento(const std::complex<T>* original_obj,const entitySet en, Memento< std::complex<T> >& memento){
      int gi=0;
      int* temp_framing =new int[en.size()] ;
      Memento< std::complex<T> >::info_type temp_info=0;
      FORALL(en,ii){
	temp_info += 2;
	temp_framing[gi]= 1;
	gi++;
      }ENDFORALL
	 memento.set_info(temp_info);
	 memento.set_framing(temp_framing);
	 Memento< std::complex<T> >::variable_type*  temp_variable= new Memento< std::complex<T> >::variable_type[2*memento.get_info()];
      int ti=0;
      gi=0;
      FORALL(en,ii){
	for(int g=0;g<temp_framing[gi];g++){
	  temp_variable[ti]=original_obj[ii].imag();
	  temp_variable[++ti]=original_obj[ii].real();
	  ti++;
	}
	gi++;
      }ENDFORALL
	 memento.set_variable(temp_variable);
  };
  
  //-----------------------set store memento-----------------------//
  //---------------------------------------------------------//
  template <class W> void set_store_memento( W* original_obj,const entitySet en, Memento<W>& memento){
    cerr<<"set_store_memento<class T> should not be called"<<endl;
  }

 template <class T> inline void set_store_memento( std::vector<T>* original_obj,const entitySet en, Memento< std::vector<T> >& memento){
   Memento< std::vector<T> >::variable_type* temp_data = memento.get_variable();
   int* temp_framing = memento.get_framing();

   //restore data in original object
   int gi=0,ti=0;
    FORALL(en,ii){
	for(int g=0;g<temp_framing[gi];g++){
	  original_obj[ii].push_back(temp_data[ti]);
	  ti++;
	}
	gi++;
    }ENDFORALL       
 };
 
 template <class T> inline void set_store_memento( std::complex<T>* original_obj,const entitySet en, Memento< std::complex<T> >& memento){
   Memento< std::complex<T> >::variable_type* temp_data = memento.get_variable();
   int* temp_framing = memento.get_framing();

   //restore data in original object
   int gi=0,ti=0;
    FORALL(en,ii){
	for(int g=0;g<temp_framing[gi];g++){
	  std::complex<T> x(temp_data[ti],temp_data[ti++]);
	  original_obj[ii] = x;
	  ti++;
	}
	gi++;
      }ENDFORALL
 };
 
 //-----------------create storeVec memento-----------------------//
 //---------------------------------------------------------------//
 template <class W> void create_storeVec_memento(const W* original_obj,const entitySet en,int size, Memento<W>& memento){
    cerr<<"create_storeVec_memento<class T> should not be called"<<endl;
  }

 template <class T> inline void create_storeVec_memento(const std::vector<T>* original_obj,const entitySet en,int size, Memento< std::vector<T> >& memento){
    std::vector<T> temp_data;
      int gi=0;
      int* temp_framing =new int[en.size()*size] ;
      Memento< std::vector<T> >::info_type temp_info=0;
      FORALL(en,ii){
	for(int j=0;j<size;j++){
	  temp_data=original_obj[ii*size+j];
	  temp_info += temp_data.size();
	  temp_framing[gi]=temp_data.size();
	  gi++;
	}
      }ENDFORALL
	 memento.set_info(temp_info);
	 memento.set_framing(temp_framing);
	 Memento< std::vector<T> >::variable_type*  temp_variable= new Memento< std::vector<T> >::variable_type[memento.get_info()];
      int ti=0;
      gi=0;
      FORALL(en,ii){
	for(int j=0;j<size;j++){
	  temp_data=original_obj[ii*size+j];
	  for(int g=0;g<temp_framing[gi];g++){
	    temp_variable[ti]=temp_data[g];
	    ti++;
	  }
	  gi++;
	}
      }ENDFORALL
	 memento.set_variable(temp_variable);
  };
 
 template <class T> inline void create_storeVec_memento(const std::complex<T>* original_obj,const entitySet en,int size, Memento< std::complex<T> >& memento){
      int gi=0;
      int* temp_framing =new int[en.size()*size] ;
      Memento< std::complex<T> >::info_type temp_info=0;
      FORALL(en,ii){
	for(int j=0;j<size;j++){
	  temp_info += 2;
	  temp_framing[gi]=1;
	  gi++;
	}
      }ENDFORALL
	 memento.set_info(temp_info);
	 memento.set_framing(temp_framing);
	 Memento< std::complex<T> >::variable_type*  temp_variable= new Memento< std::complex<T> >::variable_type[2*memento.get_info()];
      int ti=0;
      gi=0;
      FORALL(en,ii){
	for(int j=0;j<size;j++){
	  for(int g=0;g<temp_framing[gi];g++){
	    temp_variable[ti]=original_obj[ii*size+j].imag();
	    temp_variable[++ti]=original_obj[ii*size+j].real();
	    ti++;
	  }
	  gi++;
	}
      }ENDFORALL
	 memento.set_variable(temp_variable);
  };
 
  //-----------------------set storeVec memento-----------------------//
  //---------------------------------------------------------//
  template <class W> void set_storeVec_memento( W* original_obj,const entitySet en,int size, Memento<W>& memento){
    cerr<<"set_storeVec_memento<class T> should not be called"<<endl;
  };

  template <class T> inline void set_storeVec_memento( std::vector<T>* original_obj,const entitySet en,int size, Memento< std::vector<T> >& memento){
    Memento< std::vector<T> >::variable_type* temp_data = memento.get_variable();
   int* temp_framing = memento.get_framing();

   //restore data in original object
   int gi=0,ti=0;
   FORALL(en,ii){
     for(int j=0;j<size;j++){
       for(int g=0;g<temp_framing[gi];g++){
	 original_obj[ii*size+j].push_back(temp_data[ti]);
	 ti++;
       }
       gi++;
     }
   }ENDFORALL
  };
  
  template <class T> inline void set_storeVec_memento( std::complex<T>* original_obj,const entitySet en,int size, Memento< std::complex<T> >& memento){
    Memento< std::complex<T> >::variable_type* temp_data = memento.get_variable();
   int* temp_framing = memento.get_framing();

   //restore data in original object
   int gi=0,ti=0;
   FORALL(en,ii){
     for(int j=0;j<size;j++){
       for(int g=0;g<temp_framing[gi];g++){
	 std::complex<T> x(temp_data[ti],temp_data[ti++]);
	 original_obj[ii*size+j]= x;
	 ti++;
       }
       gi++;
     }
   }ENDFORALL
  };
  
//-------------------create param memento-------------------------//
  //----------------------------------------------------------//
  template <class W> void create_param_memento(const W& original_obj,const entitySet en, Memento<W>& memento){
    cerr<<"create_param_memento<class T> should not be called"<<endl;
  };

 template <class T> inline void create_param_memento(const std::vector<T>& original_obj,const entitySet en, Memento< std::vector<T> >& memento){
    std::vector<T> temp_data;
    temp_data=original_obj;

    Memento< std::vector<T> >::info_type temp_info=temp_data.size();
    memento.set_info(temp_info);

    //create dummy framing
    int* temp_framing = new int[en.size()];

    Memento< std::vector<T> >::variable_type*  temp_variable= new Memento< std::vector<T> >::variable_type[memento.get_info()];
    
    for(int g=0;g<temp_info;g++){
      temp_variable[g]=temp_data[g];
    }  
    for(int i=0;i<en.size();i++){
      temp_framing[i]=temp_info;
    }
    memento.set_variable(temp_variable);
    memento.set_framing(temp_framing);
 };
 
 template <class T> inline void create_param_memento(const std::complex<T>& original_obj,const entitySet en, Memento< std::complex<T> >& memento){
    Memento< std::complex<T> >::info_type temp_info=1;
    memento.set_info(temp_info);

    //create dummy framing
    int* temp_framing = new int[en.size()];

    Memento< std::complex<T> >::variable_type*  temp_variable= new Memento< std::complex<T> >::variable_type[2*memento.get_info()];
    int gi=0;
    for(int g=0;g<temp_info;g++){
      temp_variable[gi]=original_obj.imag();
      temp_variable[++gi]=original_obj.real();
      gi++;
    }  
    for(int i=0;i<en.size();i++){
      temp_framing[i]=temp_info;
    }
    memento.set_variable(temp_variable);
    memento.set_framing(temp_framing);
 };
 
 //--------------------set param memento----------------//
 //-----------------------------------------------------//
 template <class W> void set_param_memento( W& original_obj,const entitySet en, Memento<W>& memento){
   cerr<<"set_param_memento<class T> should not be called"<<endl;
 };

 template <class T> inline void set_param_memento( std::vector<T>& original_obj,const entitySet en, Memento< std::vector<T> >& memento){
  Memento< std::vector<T> >::variable_type* temp_data = memento.get_variable();
   int* temp_framing = memento.get_framing();

   //restore data in original object
   int gi=0;
   for(int g=0;g<temp_framing[gi];g++){
     original_obj.push_back(temp_data[g]);
   }
 };
 
 template <class T> inline void set_param_memento( std::complex<T>& original_obj,const entitySet en, Memento< std::complex<T> >& memento){
   Memento< std::complex<T> >::variable_type* temp_data = memento.get_variable();
   int* temp_framing = memento.get_framing();

   //restore data in original object
   int gi=0;
   for(int g=0;g<temp_framing[gi];g++){
     std::complex<T> x(temp_data[gi],temp_data[gi++]);
     original_obj = x;
     gi++;
   }
 };
 
 //--------------------memento hdf5write----------------//
 //-----------------------------------------------------//
 template <class T, class W> void memento_hdf5write(H5::Group group, T t, const Memento<W> &memento,const entitySet &en, int size){
   cerr<<" memento_hdf5write:This should not be called!"<<endl;
 }; 

 template <class W> void memento_hdf5write(H5::Group group, Need_Variable_data t, const Memento<W> &memento,const entitySet &en,int size){
 typedef hdf5_schema_converter_traits <W> _schema_traits_type;
   hsize_t dimf_POD[1],dimf_variable[1],dimf_framing[1];
   int RANK=1;
   dimf_POD[0]= 1;
   Memento<W>::info_type temp_info = memento.get_info();
   dimf_variable[0]=memento.get_info();
   dimf_framing[0]=en.size()*size;

    try{
      //write data
      H5::DataType variable_datatype=_schema_traits_type::get_variable_HDF5_type();
      H5::DataSpace variable_dataspace( RANK, dimf_variable );
      H5::DataSet variable_dataset = group.createDataSet( "variable", variable_datatype, variable_dataspace);
      variable_dataset.write( memento.get_variable(), variable_datatype );

      H5::DataType POD_datatype=_schema_traits_type::get_POD_HDF5_type();
      H5::DataSpace POD_dataspace( RANK, dimf_POD );
      H5::DataSet POD_dataset = group.createDataSet( "POD", POD_datatype, POD_dataspace);
      POD_dataset.write( &temp_info, POD_datatype );

      H5::DataType framing_datatype=H5::PredType::NATIVE_INT;
      H5::DataSpace framing_dataspace( RANK, dimf_framing );
      H5::DataSet framing_dataset = group.createDataSet( "framing", framing_datatype, framing_dataspace);
      framing_dataset.write( memento.get_framing(), framing_datatype );
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
 } 

 template <class W> void memento_hdf5write(H5::Group group, No_Need_Variable_data t, const Memento<W> &memento,const entitySet &en,int size){
   typedef hdf5_schema_converter_traits <W> _schema_traits_type;
   hsize_t dimf_POD[1];
   int RANK=1;
   dimf_POD[0]= 1;
   Memento<W>::info_type temp_info = memento.get_info();
    try{
      //write data,only POD data is needed to write out    
      H5::DataType POD_datatype=_schema_traits_type::get_POD_HDF5_type();
      H5::DataSpace POD_dataspace( RANK, dimf_POD );
      H5::DataSet POD_dataset = group.createDataSet( "POD", POD_datatype, POD_dataspace);
      POD_dataset.write( &temp_info, POD_datatype );
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
 };

  //-----------------------memento hdf5read--------------------//
 //------------------------------------------------------------//
 template <class W> void memento_hdf5read(H5::Group group, Need_Variable_data t, Memento<W> &memento,entitySet &en){
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
	Memento<W>::variable_type *data_variable = new Memento<W>::variable_type[dims_variable[0]];
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
};

  template <class W> void memento_hdf5read(H5::Group group, No_Need_Variable_data t, Memento<W> &memento,entitySet &en){
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
  };

}
#endif
