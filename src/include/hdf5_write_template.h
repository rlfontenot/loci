#ifndef HDF5_WRITE_TEMPLATE_H_
#define HDF5_WRITE_TEMPLATE_H_

#include <hdf5_memento.h>
#include <hdf5_traits.h>
#include <hdf5CC/H5cpp.h>
#include <Tools/stream.h>
#include <typeinfo>
namespace Loci {

 //------------------domain hdf5write--------------------//
 //------------------------------------------------------//
  inline void domain_hdf5write(H5::Group group, const entitySet& en){
    hsize_t dimf_domain[1];
    int RANK=1; 
    int num_intervals=en.num_intervals();
    dimf_domain[0]=num_intervals*2;
    interval *it = new interval[num_intervals];
    int *data_domain = new int[num_intervals*2];//get the domain data
    for(int i=0;i<num_intervals;i++){
      it[i]=en[i];
      data_domain[i*2]=it[i].first;
      data_domain[i*2+1]=it[i].second;
    }
    try{
      //write domain
      H5::DataSpace dataspace_domain( RANK, dimf_domain );
      H5::DataSet dataset_domain = group.createDataSet( "domain", H5::PredType::NATIVE_INT, dataspace_domain );
      dataset_domain.write( data_domain, H5::PredType::NATIVE_INT );
      
      delete [] it;
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
  };

  //------------------------store hdf5write--------------------//
  //-----------------------------------------------------------//
  template <class T,class W>  void store_hdf5write(H5::Group group,const W* base_ptr,const entitySet &en) {
    cerr<<"store_hdf5write:This should not be called!"<<endl;
  };

 template <class W> void store_hdf5write(H5::Group group,DEFAULT_CONVERTER g,const W* base_ptr,const entitySet &en) {
   
   int RANK=1;
   hsize_t dimf_store[1];
   std::ostringstream oss;
   //cout<<"store write out "<<endl;
   oss << '{' << en << std::endl ;
   FORALL(en,ii) {
     oss << base_ptr[ii] << std::endl ;
     //cout<< base_ptr[ii] <<endl;
   }ENDFORALL ;
   oss << '}' << std::endl ;
   
   std::string memento=oss.str();
   hsize_t size=memento.length();
   dimf_store[0]= size+1;
   try{
      H5::DataSpace dataspace( RANK, dimf_store );
      H5::DataSet dataset = group.createDataSet( "store",H5::PredType::NATIVE_CHAR , dataspace);
      dataset.write( memento.c_str(), H5::PredType::NATIVE_CHAR );
   }
   catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
   catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
   catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
   
 }

  template <class W> void store_hdf5write(H5::Group group,IDENTITY_CONVERTER g,const W* base_ptr,const entitySet &en) {
    
      typedef hdf5_schema_traits<W> _schema_traits_type;
      int dim=en.size();
      hsize_t dimf[1];
      dimf[0]=dim;
      int RANK=1;

      int num_intervals=en.num_intervals();
      interval *it = new interval[num_intervals];
      int bound=0;
      for(int i=0;i<num_intervals;i++){
        it[i]=en[i];
      }

      if(en.Min()<0&&en.Max()>0){
	bound=en.Max()-en.Min()+1;
      }
      else if(en.Min()<0){
	bound=abs(en.Min());
      }
      else
	bound=en.Max();
	//for negative domain, we do the coord. transformation
      if(en.Min()<0){
	for(int i=0;i<num_intervals;i++){
	  it[i].first=it[i].first+abs(en.Min());
	  it[i].second=it[i].second+abs(en.Min());
	}
      }

      //declare the variables used by hyperslab
      hsize_t dim_mem[1];
      dim_mem[0]=bound;
      hssize_t start_mem[1];
      hsize_t stride_mem[1];
      hsize_t count_mem[1];
      hsize_t block_mem[1];
      hssize_t start_file[1];
      start_file[0]=0;
      stride_mem[0]=1;	    
      block_mem[0]=1;

      domain_hdf5write(group,en);
      try{
	//create the store data
	H5::DataType datatype=_schema_traits_type::get_type();
	H5::DataSpace dataspace_store( RANK, dimf );
	H5::DataSet dataset_store = group.createDataSet( "store", datatype, dataspace_store );
	//get data using HDF5 hyperslab	
	H5::DataSpace dataspace_memory(RANK,dim_mem);
	for(int i=0;i<num_intervals;i++){
	    start_mem[0]=it[i].first;
	    count_mem[0]=it[i].second-it[i].first+1;
	    dataspace_memory.selectHyperslab(H5S_SELECT_SET,count_mem,start_mem,stride_mem,block_mem);	
	    dataspace_store.selectHyperslab(H5S_SELECT_SET,count_mem,start_file,stride_mem,block_mem);
	    start_file[0]=start_file[0]+count_mem[0];//for next interval
	    dataset_store.write(base_ptr,datatype,dataspace_memory,dataspace_store);
	}
      }
      catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
      catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();} 
      catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}

      //reclaim memory
      delete [] it;
  };

 template <class W> void store_hdf5write(H5::Group group,USER_DEFINED_CONVERTER g,const W* base_ptr,const entitySet &en) {   
   typedef hdf5_schema_converter_traits<W> schema_converter;
   typedef typename hdf5_schema_converter_traits<W>::need_variable_selector need_selector;

   //get memento information
   Memento<W> memento;
   create_store_memento(base_ptr,en,memento);

   //write out the domain   
    domain_hdf5write(group,en);

   //write out the data to HDF5 file
   need_selector n_selector;
   memento_hdf5write(group,n_selector, memento,en,1);
   
 };
 
 //--------------get store domain------------------------//
  //------------------------------------------------------//
  template <class T> int get_storeVec_size(H5::Group group,T g){
    cerr<<" get_storeVec_size:This should not be called!"<<endl;
  };

  template <> inline int get_storeVec_size(H5::Group group,DEFAULT_CONVERTER g){
    char ch;
    int size;
   try{
     H5::DataSet dataset_store = group.openDataSet( "store");
     H5::DataSpace dataspace_store = dataset_store.getSpace();
     hsize_t dims_store[1];
     dataspace_store.getSimpleExtentDims( dims_store, NULL);
     char* memento= new char[dims_store[0]];
     dataset_store.read( memento, H5::PredType::NATIVE_CHAR );
     std::istringstream iss(memento);
     
     do ch = iss.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '{') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        iss.putback(ch) ;
      }
      
      entitySet e ;
      iss >> e ;
      iss >> size ;

      delete [] memento;
   }
   catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
   catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();} 
   catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}

   return size;
  };

  template <> inline int get_storeVec_size(H5::Group group,IDENTITY_CONVERTER g){
    int size;
    try{
      //read in size
      H5::DataSet dataset_size = group.openDataSet( "size");
      dataset_size.read(&size,H5::PredType::NATIVE_INT);
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    return size;
  };

 template <> inline int get_storeVec_size(H5::Group group,USER_DEFINED_CONVERTER g){
    IDENTITY_CONVERTER j;
    int num= get_storeVec_size(group,j);
    return num;
  };
  //--------------get store domain------------------------//
  //------------------------------------------------------//
  template <class T> entitySet get_store_domain(H5::Group group,T g){
     cerr<<" get_store_domain:This should not be called!"<<endl;
  };
  
  template <> inline entitySet get_store_domain(H5::Group group,DEFAULT_CONVERTER g){
    char ch;
   entitySet e;
   try{
     H5::DataSet dataset_store = group.openDataSet( "store");
     H5::DataSpace dataspace_store = dataset_store.getSpace();
     hsize_t dims_store[1];
     dataspace_store.getSimpleExtentDims( dims_store, NULL);
     char* data_store= new char[dims_store[0]];
     dataset_store.read( data_store, H5::PredType::NATIVE_CHAR );
     std::istringstream iss(data_store);
    do ch = iss.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      iss.putback(ch) ;
    }
    iss >> e ;

    delete [] data_store;
   }
   catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
   catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();} 
   catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}

   return e;
}
  
  template <> inline entitySet get_store_domain(H5::Group group,IDENTITY_CONVERTER h){
    entitySet num;
    try{
	//get domain data
	H5::DataSet dataset_domain = group.openDataSet( "domain");
	H5::DataSpace dataspace_domain = dataset_domain.getSpace();
	hsize_t dims_domain[1];
	dataspace_domain.getSimpleExtentDims( dims_domain, NULL);
	int *data_domain = new int[dims_domain[0]];
	dataset_domain.read( data_domain, H5::PredType::NATIVE_INT );
	for(int i=0;i<dims_domain[0];i++){
	  num |=interval(data_domain[i],data_domain[i+1]);
	  i++;
	}
	delete [] data_domain;
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}

    return num;
  };
  
  template <> inline entitySet get_store_domain(H5::Group group,USER_DEFINED_CONVERTER h){
    IDENTITY_CONVERTER j;
    entitySet num= get_store_domain(group,j);
    return num;
  };

  //-------------------store hdf5read----------------//
  //-------------------------------------------------//
  template <class T,class W> void store_hdf5read(H5::Group group,T t,W* base_ptr,entitySet &en){
    cerr<<" store_hdf5read:This should not be called!"<<endl;
    };

  template <class W> void store_hdf5read(H5::Group group,DEFAULT_CONVERTER g,W* base_ptr,entitySet &en){
   char ch;
   entitySet e;
   try{
     H5::DataSet dataset_store = group.openDataSet( "store");
     H5::DataSpace dataspace_store = dataset_store.getSpace();
     hsize_t dims_store[1];
     dataspace_store.getSimpleExtentDims( dims_store, NULL);
     char* memento= new char[dims_store[0]];
     dataset_store.read( memento, H5::PredType::NATIVE_CHAR );
     std::istringstream iss(memento);

     do ch = iss.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      iss.putback(ch) ;
    }
    iss >> e ;
    //cout<<"store read in "<<endl;
    FORALL(en,ii) {
      iss >> base_ptr[ii] ;
      //cout << base_ptr[ii] <<endl;
    } ENDFORALL ;
        
    do ch = iss.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      iss.putback(ch) ;
    }
   delete [] memento;
   }
   catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
   catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
   catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
 }

 template <class W> void store_hdf5read(H5::Group group,IDENTITY_CONVERTER g,W* base_ptr,entitySet &en){
   typedef hdf5_schema_traits<W> _schema_traits_type;
    hsize_t dims_store[1];
      
      try{
	//get store data  
	H5::DataSet dataset_store = group.openDataSet( "store");
	H5::DataSpace dataspace_store = dataset_store.getSpace();
	dataspace_store.getSimpleExtentDims( dims_store, NULL);
	int RANK = dataspace_store.getSimpleExtentNdims();

	//set the intervals
	int num_intervals=en.num_intervals();
	interval *it = new interval[num_intervals];
	int bound=0;
	for(int i=0;i<num_intervals;i++){
	  it[i]=en[i];
	}

	if(en.Min()<0&&en.Max()>0){
	bound=en.Max()-en.Min()+1;
      }
      else if(en.Min()<0){
	bound=abs(en.Min());
      }
      else
	bound=en.Max();
	//for negative domain, we do the coord. transformation
      if(en.Min()<0){
	for(int i=0;i<num_intervals;i++){
	  it[i].first=it[i].first+abs(en.Min());
	  it[i].second=it[i].second+abs(en.Min());
	}
      }

	//declear the variables used by hyperslab
	hsize_t dim_mem[1];
	dim_mem[0]=bound;
	hssize_t start_mem[1];
	hsize_t stride_mem[1];
	hsize_t count_mem[1];
	hsize_t block_mem[1];
	hssize_t start_file[1];
	start_file[0]=0;
	stride_mem[0]=1;	    
	block_mem[0]=1;

	//get data using HDF5 hyperslab	
	H5::DataSpace dataspace_memory(RANK,dim_mem);
	H5::DataType datatype=_schema_traits_type::get_type();
	for(int i=0;i<num_intervals;i++){
	    start_mem[0]=it[i].first;
	    count_mem[0]=it[i].second-it[i].first+1;
	    dataspace_memory.selectHyperslab(H5S_SELECT_SET,count_mem,start_mem,stride_mem,block_mem);	
	    dataspace_store.selectHyperslab(H5S_SELECT_SET,count_mem,start_file,stride_mem,block_mem);
	    start_file[0]=start_file[0]+count_mem[0];//for next interval
	    dataset_store.read(base_ptr,datatype,dataspace_memory,dataspace_store);
	}

      //reclaim memory
      delete [] it;
      } 
      catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
      catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
      catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
 }
  
 template <class W> void store_hdf5read(H5::Group group,USER_DEFINED_CONVERTER g, W* base_ptr, entitySet &en) {
   typedef hdf5_schema_converter_traits<W> schema_converter;
   typedef typename hdf5_schema_converter_traits<W>::need_variable_selector need_selector;

   //get memento information
   Memento<W> memento;
   //domain data already get from get_store_domain()

   //get data
   need_selector n_selector;
   memento_hdf5read(group,n_selector,memento,en);
   set_store_memento(base_ptr,en,memento);
 };


  //------------------------storeVec hdf5write--------------------//
  //-----------------------------------------------------------//
  template <class T,class W>  void storeVec_hdf5write(H5::Group group,const W* base_ptr,const entitySet &en,const int size) {
    cerr<<"storeVec_hdf5write:This should not be called!"<<endl;
  };

 template <class W>  void storeVec_hdf5write(H5::Group group,DEFAULT_CONVERTER g, W* base_ptr,const entitySet &en,const int size) {
   int RANK=1;
   hsize_t dimf_store[1];
   std::ostringstream oss;

   oss << '{' << en << std::endl ;
   oss << size << std::endl ;
   //cout<<"storeVec write out"<<endl;
   FORALL(en,ii) {
     W* p = base_ptr + ii*size ;
     for(int i=0;i<size;++i,++p)
       oss << *p << " " ;
     oss << std::endl ;
     //cout<<*p<<endl;
   }ENDFORALL ;
   oss << '}' << std::endl ;
   
   std::string memento=oss.str();
   hsize_t m_size=memento.length();
   dimf_store[0]= m_size+1;
   try{
      H5::DataSpace dataspace( RANK, dimf_store );
      H5::DataSet dataset = group.createDataSet( "store",H5::PredType::NATIVE_CHAR , dataspace);
      dataset.write( memento.c_str(), H5::PredType::NATIVE_CHAR );
   }
   catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
   catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
   catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
  };

 template <class W>  void storeVec_hdf5write(H5::Group group,IDENTITY_CONVERTER g,const W* base_ptr,const entitySet &en,const int size) {
   typedef hdf5_schema_traits<W> _schema_traits_type;
      int dim=en.size();
      hsize_t dimf[1];
      dimf[0]=dim*size;
      int RANK=1;

      int num_intervals=en.num_intervals();
      interval *it = new interval[num_intervals];
      int bound=0;
      for(int i=0;i<num_intervals;i++){
        it[i]=en[i];
      }

      if(en.Min()<0&&en.Max()>0){
	bound=en.Max()-en.Min()+1;
      }
      else if(en.Min()<0){
	bound=abs(en.Min());
      }
      else
	bound=en.Max();
	//for negative domain, we do the coord. transformation
      if(en.Min()<0){
	for(int i=0;i<num_intervals;i++){
	  it[i].first=it[i].first+abs(en.Min());
	  it[i].second=it[i].second+abs(en.Min());
	}
      }

      //declare the variables used by hyperslab
      hsize_t dim_mem[1];
      dim_mem[0]=(bound+1)*size;
      hssize_t start_mem[1];
      hsize_t stride_mem[1];
      hsize_t count_mem[1];
      hsize_t block_mem[1];
      hssize_t start_file[1];
      start_file[0]=0;
      stride_mem[0]=1;	    
      block_mem[0]=1;
      hsize_t dim_size[1]={1};

      domain_hdf5write(group,en);
      try{
	//write out size
	H5::DataSpace dataspace_size( RANK, dim_size );
	H5::DataSet dataset_size = group.createDataSet( "size", H5::PredType::NATIVE_INT, dataspace_size );
	 dataset_size.write(&size,H5::PredType::NATIVE_INT);
	//create the store data
	H5::DataType datatype=_schema_traits_type::get_type();
	H5::DataSpace dataspace_store( RANK, dimf );
	H5::DataSet dataset_store = group.createDataSet( "store", datatype, dataspace_store );
	//get data using HDF5 hyperslab	
	H5::DataSpace dataspace_memory(RANK,dim_mem);
	for(int i=0;i<num_intervals;i++){
	  start_mem[0]=(it[i].first)*size;
	  count_mem[0]=(it[i].second-it[i].first+1)*size;
	  dataspace_memory.selectHyperslab(H5S_SELECT_SET,count_mem,start_mem,stride_mem,block_mem);	
	  dataspace_store.selectHyperslab(H5S_SELECT_SET,count_mem,start_file,stride_mem,block_mem);
	  start_file[0]=start_file[0]+count_mem[0];//for next interval
	  dataset_store.write(base_ptr,datatype,dataspace_memory,dataspace_store);
	}
      }
      catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
      catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();} 
      catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}

      //reclaim memory
      delete [] it;
  };

 template <class W>  void storeVec_hdf5write(H5::Group group,USER_DEFINED_CONVERTER g,const W* base_ptr,const entitySet &en,const int size) {
    typedef hdf5_schema_converter_traits<W> schema_converter;
   typedef typename hdf5_schema_converter_traits<W>::need_variable_selector need_selector;

   hsize_t dim_size[1]={1};
   int RANK=1;
   try{
     //write out size
     H5::DataSpace dataspace_size( RANK, dim_size );
     H5::DataSet dataset_size = group.createDataSet( "size", H5::PredType::NATIVE_INT, dataspace_size );
     dataset_size.write(&size,H5::PredType::NATIVE_INT);
   }
   catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
   catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();} 
   catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}


   //get memento information
   Memento<W> memento;
   create_storeVec_memento(base_ptr,en,size,memento);

   //write out the domain   
    domain_hdf5write(group,en);

   //write out the data to HDF5 file
   need_selector n_selector;
   memento_hdf5write(group,n_selector, memento,en,size);
   
  };

 //-------------------storeVec hdf5read----------------//
  //-------------------------------------------------//
  template <class T,class W> void storeVec_hdf5read(H5::Group group,T t,W* base_ptr,entitySet &en,int size){
    cerr<<" store_hdf5read:This should not be called!"<<endl;
  };

 template <class W> void storeVec_hdf5read(H5::Group group,DEFAULT_CONVERTER t,W* base_ptr,entitySet &en,int size){
   char ch;
   try{
     H5::DataSet dataset_store = group.openDataSet( "store");
     H5::DataSpace dataspace_store = dataset_store.getSpace();
     hsize_t dims_store[1];
     dataspace_store.getSimpleExtentDims( dims_store, NULL);
     char* memento= new char[dims_store[0]];
     dataset_store.read( memento, H5::PredType::NATIVE_CHAR );
     std::istringstream iss(memento);
     
     do ch = iss.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '{') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        iss.putback(ch) ;
      }
      
      entitySet e ;
      iss >> e ;
      iss >> size ;
      //cout<<"storeVec read in"<<endl;
      FORALL(e,ii) {
        W * p = base_ptr + ii*size ;
        for(int i=0;i<size;++i,++p)
          iss >> *p;
	//cout<<*p<<endl;
      } ENDFORALL ;
      
      do ch = iss.get(); while(ch==' ' || ch=='\n') ;
      if(ch != '}') {
        std::cerr << "Incorrect Format while reading store" << std::endl ;
        iss.putback(ch) ;
	}
   }
   catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
   catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
   catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
  };

 template <class W> void storeVec_hdf5read(H5::Group group,IDENTITY_CONVERTER t,W* base_ptr,entitySet &en,int size){
    typedef hdf5_schema_traits<W> _schema_traits_type;
    hsize_t dims_store[1];
    
    size=get_storeVec_size(group,t);
      try{
	//get store data  
	H5::DataSet dataset_store = group.openDataSet( "store");
	H5::DataSpace dataspace_store = dataset_store.getSpace();
	dataspace_store.getSimpleExtentDims( dims_store, NULL);
	int RANK = dataspace_store.getSimpleExtentNdims();

	//set the intervals
	int num_intervals=en.num_intervals();
	interval *it = new interval[num_intervals];
	int bound=0;
	for(int i=0;i<num_intervals;i++){
	  it[i]=en[i];
	}

	if(en.Min()<0&&en.Max()>0){
	  bound=en.Max()-en.Min()+1;
	}
	else if(en.Min()<0){
	  bound=abs(en.Min());
	}
	else
	  bound=en.Max();
	//for negative domain, we do the coord. transformation
	if(en.Min()<0){
	  for(int i=0;i<num_intervals;i++){
	    it[i].first=it[i].first+abs(en.Min());
	    it[i].second=it[i].second+abs(en.Min());
	  }
	}
	
	//declear the variables used by hyperslab
	hsize_t dim_mem[1];
	dim_mem[0]=(bound+1)*size;
	hssize_t start_mem[1];
	hsize_t stride_mem[1];
	hsize_t count_mem[1];
	hsize_t block_mem[1];
	hssize_t start_file[1];
	start_file[0]=0;
	stride_mem[0]=1;	    
	block_mem[0]=1;

	//get data using HDF5 hyperslab	
	H5::DataSpace dataspace_memory(RANK,dim_mem);
	H5::DataType datatype=_schema_traits_type::get_type();
	for(int i=0;i<num_intervals;i++){
	    start_mem[0]=(it[i].first)*size;
	    count_mem[0]=(it[i].second-it[i].first+1)*size;
	    dataspace_memory.selectHyperslab(H5S_SELECT_SET,count_mem,start_mem,stride_mem,block_mem);	
	    dataspace_store.selectHyperslab(H5S_SELECT_SET,count_mem,start_file,stride_mem,block_mem);
	    start_file[0]=start_file[0]+count_mem[0];//for next interval
	    dataset_store.read(base_ptr,datatype,dataspace_memory,dataspace_store);
	}

      //reclaim memory
      delete [] it;
      } 
      catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
      catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
      catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
 };

 template <class W> void storeVec_hdf5read(H5::Group group,USER_DEFINED_CONVERTER t,W* base_ptr,entitySet &en,int size){
   typedef hdf5_schema_converter_traits<W> schema_converter;
   typedef typename hdf5_schema_converter_traits<W>::need_variable_selector need_selector;

   //get memento information
   Memento<W> memento;
   //domain data already get from get_store_domain()

   //get data
   need_selector n_selector;
   memento_hdf5read(group,n_selector,memento,en);
   set_storeVec_memento(base_ptr,en,size,memento);
  };

 //------------------parameter hdf5write-----------------------//
 //-----------------------------------------------------------//
  template <class T,class W>  void param_hdf5write(H5::Group group,const W &param,const entitySet &en) {
    cerr<<"hdf5write:This should not be called!"<<endl;
  }; 

  template <class W> void param_hdf5write(H5::Group group,DEFAULT_CONVERTER g,const W &param,const entitySet &en) {
    typedef hdf5_schema_traits <W> _schema_traits_type;
    std::ostringstream oss;
    oss<<param;
    //cout<<"param write out "<<endl;
    //cout<<param<<endl;
    std::string memento=oss.str();
    hsize_t size=memento.length();

    hsize_t dimf_param[1];
    int RANK=1;    
    int num_intervals=en.num_intervals();
    dimf_param[0]= size+1;
    interval *it = new interval[num_intervals];

    domain_hdf5write(group,en);
    try{
      //write parameter
      H5::DataType datatype=_schema_traits_type::get_type();
      H5::DataSpace dataspace( RANK, dimf_param );
      H5::DataSet dataset = group.createDataSet( "param", datatype, dataspace);
      dataset.write( memento.c_str(), datatype );
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}

    delete [] it;
  };

  template <class W> void param_hdf5write(H5::Group group,IDENTITY_CONVERTER g,const W &param,const entitySet &en) {
   typedef hdf5_schema_traits <W> _schema_traits_type;
   hsize_t dimf_param[1];
   int RANK=1;
   int num_intervals=en.num_intervals();
   dimf_param[0]= 1;
   interval *it = new interval[num_intervals];

   domain_hdf5write(group,en);
    try{
      //write parameter
      H5::DataType datatype=_schema_traits_type::get_type();
      H5::DataSpace dataspace( RANK, dimf_param );
      H5::DataSet dataset = group.createDataSet( "param", datatype, dataspace);
      dataset.write( &param, datatype );
    }
    catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
    catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
    catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}

    delete [] it;
  };

  template <class W> void param_hdf5write(H5::Group group,USER_DEFINED_CONVERTER g,const W &param,const entitySet &en) {
   typedef hdf5_schema_converter_traits<W> schema_converter;
   typedef typename hdf5_schema_converter_traits<W>::need_variable_selector need_selector;

   //get memento information
   Memento<W> memento;
   create_param_memento(param,en,memento);

   //write out the domain   
    domain_hdf5write(group,en);

   //write out the data to HDF5 file
   need_selector n_selector;
   memento_hdf5write(group,n_selector, memento,en,1);
  };

  //---------------------parameter hdf5read---------------------//
  //------------------------------------------------------------//
 template <class T,class W> void param_hdf5read(H5::Group group,T t,W &param){
    cerr<<" hdf5read:This should not be called!"<<endl;
    };

 template <class W> entitySet param_hdf5read(H5::Group group,DEFAULT_CONVERTER g,W &param){
   typedef hdf5_schema_traits <W> _schema_traits_type;
   entitySet num;	

   try{
     //get domain data
     H5::DataSet dataset_domain = group.openDataSet( "domain");
     H5::DataSpace dataspace_domain = dataset_domain.getSpace();
     hsize_t dims_domain[1];
     dataspace_domain.getSimpleExtentDims( dims_domain, NULL);
     int *data_domain = new int[dims_domain[0]];
     dataset_domain.read( data_domain, H5::PredType::NATIVE_INT );
     for(int i=0;i<dims_domain[0];i++){
       num |=interval(data_domain[i],data_domain[i+1]);
       i++;
     }

     //get param data
     H5::DataType datatype=_schema_traits_type::get_type();
     H5::DataSet dataset_param = group.openDataSet( "param");
     H5::DataSpace dataspace_param = dataset_param.getSpace();
     hsize_t dims_param[1];
     dataspace_param.getSimpleExtentDims( dims_param, NULL);
     char* data_param= new char[dims_param[0]];
     dataset_param.read( data_param, datatype );
     std::istringstream iss(data_param);
     iss >> param;
     //cout<<"param read in "<<endl;
     //cout<<param<<endl;
     delete [] data_param;
     delete [] data_domain;
   }
   catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
   catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
   catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
   
   return num;
 };

 template <class W> entitySet param_hdf5read(H5::Group group,IDENTITY_CONVERTER g,W &param){ 
   typedef hdf5_schema_traits <W> _schema_traits_type;
   entitySet num;
   try{
     //get domain data
     H5::DataSet dataset_domain = group.openDataSet( "domain");
     H5::DataSpace dataspace_domain = dataset_domain.getSpace();
     hsize_t dims_domain[1];
     dataspace_domain.getSimpleExtentDims( dims_domain, NULL);
     int *data_domain = new int[dims_domain[0]];
     dataset_domain.read( data_domain, H5::PredType::NATIVE_INT );
     for(int i=0;i<dims_domain[0];i++){
       num |=interval(data_domain[i],data_domain[i+1]);
       i++;
     }

     //get param data
     H5::DataType datatype=_schema_traits_type::get_type();
     H5::DataSet dataset_param = group.openDataSet( "param");
     H5::DataSpace dataspace_param = dataset_param.getSpace();
     hsize_t dims_param[1];
     dataspace_param.getSimpleExtentDims( dims_param, NULL);
     dataset_param.read( &param, datatype );

     delete [] data_domain;
   }
   catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
   catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
   catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}
   
   return num;
 };

 template <class W> entitySet param_hdf5read(H5::Group group,USER_DEFINED_CONVERTER g,W &param){ 
   typedef hdf5_schema_converter_traits<W> schema_converter;
   typedef typename hdf5_schema_converter_traits<W>::need_variable_selector need_selector;
   entitySet en;

   try{
     //get domain data
     H5::DataSet dataset_domain = group.openDataSet( "domain");
     H5::DataSpace dataspace_domain = dataset_domain.getSpace();
     hsize_t dims_domain[1];
     dataspace_domain.getSimpleExtentDims( dims_domain, NULL);
     int *data_domain = new int[dims_domain[0]];
     dataset_domain.read( data_domain, H5::PredType::NATIVE_INT );
     for(int i=0;i<dims_domain[0];i++){
       en |=interval(data_domain[i],data_domain[i+1]);
       i++;
     }
   }
   catch( H5::HDF5DatasetInterfaceException error ){error.printerror();}
   catch( H5::HDF5DataspaceInterfaceException error ){error.printerror();}
   catch( H5::HDF5DatatypeInterfaceException error ){error.printerror();}

   //get data
   Memento<W> memento;
   need_selector n_selector;
   memento_hdf5read(group,n_selector,memento,en);
   set_param_memento(param,en,memento);

   return en;
 }

}
#endif
