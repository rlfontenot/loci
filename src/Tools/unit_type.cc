#include <Tools/unit_type.h>

using namespace std ;

namespace Loci {
  //three tables of unit type - basic, composite, reference types----//
  UNIT_type::basic_units UNIT_type::basic_unit_table[]={
    "metre",Length,1,
    "kilogram",Mass,1,
    "second",Time,1,
    "kelvin",Temperature,1,
    "ampere",Electric_current,1,
    "mole",Amount_of_substance,1,
    "candela",Luminous_intensity,1,
  };
  UNIT_type::basic_units UNIT_type::cgs_basic_unit_table[]={
    "centimetre",Length,1,
    "gram",Mass,1,
    "second",Time,1,
    "kelvin",Temperature,1,
    "ampere",Electric_current,1,
    "mole",Amount_of_substance,1,
    "candela",Luminous_intensity,1,
    };

  UNIT_type::composite_units UNIT_type::composite_unit_table[]={
    //abbreviation of SI
    "m","metre",1,
    "kg","kilogram",1,
    "g","kilogram",0.001,
    "gram","kilogram",0.001,
    "s","second",1,
    "K","kelvin",1,
    "A","ampere",1,
    "mol","mole",1,
    "cd","candela",1,

    //metric system
    "centimetre","metre",0.01,//length
    "cm","metre",0.01,
    "kilometre","metre",1000,
    "km","metre",1000,
    "milimetre","metre",0.001,
    "mm","metre",0.001,

    "minute","second",60,//time
    "min","second",10,
    "hour","second",3600,
    "h","second",3600,
    "day","second",86400,
    "d","second",86400,
    "year","second",31536000,
    "y","second",31536000,
    "shake","second",1e-8,

    "Fahrenheit","kelvin",0.5555556,//temperature interval
    "F","kelvin",0.5555556,
    "Rankine","kelvin",0.5555556,
    "R","kelvin",0.5555556,
    "Celsius","kelvin",1,
    "C","kelvin",1,

    "m^2","m^2",1,//area

    //composite unit
    "newton","kilogram*metre/second/second",1,
    "N","kilogram*metre/second/second",1,
    "joule","kilogram*metre/second/second*metre",1,
    "J","kilogram*metre/second/second*metre",1,
    "watt","kilogram*metre/second/second*metre/second",1,
    "W","kilogram*metre/second/second*metre/second",1,
    "pascal","kilogram*metre/second/second/metre/metre",1,
    "Pa","kilogram*metre/second/second/metre/metre",1
  };

  UNIT_type::reference_units UNIT_type::reference_unit_table[]={
    "acre","m^2",404.6873,//area
    "are","m^2",100,
    "a","m^2",100,
    "barn","m^2",1e-28,
    "b","m^2",1e-28,
    "hectare","m^2",1.0e4,
    "ha","m^2",1.0e4,

    "calorie","joule",4.19002,//energy or work
    "cal","joule",4.19002,
    "electronvolt","joule",1.602177e-19,
    "eV","joule",1.602177e-19,
    "erg","joule",1.0e-7,
    "kilocalorie","joule",4.19002e3,
    "kcal","joule",4.19002e3,
    "kW*h","joule",3.6e6,
    "quad","joule",1.055056e18,

    "dyne","N",1.0e-5,//force
    "dyn","N",1.0e-5,
    "kgf","N",9.80665,
    "kp","N",9.80665,
    "kip","N",4.448222e3,
    "ozf","N",2.780139e-1,
    "poundal","N",1.38255e-1,
    "lbf","N",4.448222,
    "tonforce","N",8.896443e3,

    "chain","m",2.011684e1,//length
    "ch","m",2.011684e1,
    "fathom","m",1.828804,
    "fermi","m",1e-15,
    "feet","m",0.3048,
    "foot","m",0.3048,
    "ft","m",0.3048,
    "inch","m",2.54e-2,
    "in","m",2.54e-2,
    "microinch","m",2.54e-8,
    "micro","m",1e-6,
    "mil","m",2.54e-5,
    "mile","m",1609.344,
    "mi","m",1609.344,
    "rod","m",5.029210,
    "rd","m",5.029210,
    "yard","m",0.9144,
    "yd","m",0.9144,

    "carat","kg",2.0e-4,//mass
    "grain","kg",6.479891e-5,
    "gr","kg",6.479891e-5,
    "ounce","kg",2.834952e-2,
    "oz","kg",2.834952e-2,
    "pennyweight","kg",1.555174e-3,
    "dwt","kg",1.555174e-3,
    "pound","kg",0.4535924,
    "lb","kg",0.4535924,
    "slug","kg",1.459390e1,
    "ton","kg",1e3,

    "denier","kg/m",1.111111e-7,//mass divided by length
    "tex","kg/m",1e-6,

    "darcy","m^2",9.869233e-13,//permeability
    "perm","kg/(Pa*s*(m^2))",5.72135e-11,//(0C)

    "horsepower","W",7.354988e2,//power(metric)
    "hp","W",7.354988e2,

    "atmosphere","Pa",1.01325e5,//pressure or stress
    "bar","Pa",1e5,
    "cmHg","Pa",1.333224e3,
    "cmH2O","Pa",9.80665e1,
    "ftHg","Pa",4.063666e4,
    "ftH2O","Pa",2.989067e3,
    "inHg","Pa",3.686389e3,
    "inH2O","Pa",2.490889e2,
    "ksi","Pa",6.894757e6,
    "millibar","Pa",1e2,
    "mbar","Pa",1e2,
    "psi","Pa",6.894757e3,
    "torr","Pa",1.333224e2,
    "Torr","Pa",1.333224e2,

    "centipoise","Pa*s",1.0e-3,//viscosity,dynamic
    "cP","Pa*s",1e-3,
    "poise","Pa*s",1e-1,
    "P","Pa*s",1e-1,
    "rhe","1/(Pa*s)",1e1,
    
    "centistokes","(m^2)/s",1e-6,//viscosity,kinematic
    "stokes","(m^2)/s",1e-4,
    "St","(m^2)/s",1e-4,
    
    "acrefoot","m^3",1.233489e3,//volume
    "barrel","m^3",1.589873e-1,
    "bll","m^3",1.589873e-1,
    "bushel","m^3",3.523907e-2,
    "cord","m^3",3.624556,
    "cup","m^3",2.365882e-4,
    "gallon","m^3",4.54609e-3,
    "gal","m^3",4.54609e-3,
    "gill","m^3",1.420653e-4,
    "gi","m^3",1.420653e-4,
    "liter","m^3",1e-3,
    "peck","m^3",8.809768e-3,
    "pk","m^3",8.809768e-3,
    "pint","m^3",4.731765e-4,
    "quart","m^3",9.463529e-4,
    "stere","m^3",1,
    "st","m^3",1,
    "tablespoon","m^3",1.478676e-5,
    "teaspoon","m^3",4.928922e-6,

    "clo","(m^2)*K/W",1.55e-1,//thermal insulance

    "atm","pascal",101325,
    "bar","pascal",100000,
    "kPa","pascal",1000,
    "Btu","joule",1055.87
  };

  UNIT_type::composite_units UNIT_type::cgs_composite_unit_table[]={
    //abbreviation of SI
    "m","centimetre",0.1,
    "kilogram","gram",1000,
    "g","gram",1,
    "kg","gram",1000,
    "s","second",1,
    "K","kelvin",1,
    "A","ampere",1,
    "mol","mole",1,
    "cd","candela",1,

    //metric system
    "kilometre","centimetre",1000000,
    "km","centimetre",1000000,
    "metre","centimetre",1000,
    "milimetre","centimetre",0.1,
    "mm","centimetre",0.1,
    "cm","centimetre",1,

    "minute","second",60,//time
    "min","second",10,
    "hour","second",3600,
    "h","second",3600,
    "day","second",86400,
    "d","second",86400,
    "year","second",31536000,
    "y","second",31536000,
    "shake","second",1e-8,

    "Fahrenheit","kelvin",0.5555556,//temperature interval
    "F","kelvin",0.5555556,
    "Rankine","kelvin",0.5555556,
    "R","kelvin",0.5555556,
    "Celsius","kelvin",1,
    "C","kelvin",1,

    //composite unit
    "newton","gram*centimetre/second/second",10000,
    "N","gram*centimetre/second/second",10000,
    "joule","gram*centimetre/second/second*centimetre",100000,
    "J","gram*centimetre/second/second*centimetre",100000,
    "watt","gram*centimetre/second/second*centimetre/second",100000,
    "W","gram*centimetre/second/second*centimetre/second",100000,
    "pascal","gram*centimetre/second/second/centimetre/centimetre",10000,
    "Pa","gram*centimetre/second/second/centimetre/centimetre",10000
  };
  
  UNIT_type::default_units UNIT_type::default_unit_table[]={
    "general"," ",
    "length","metre",
    "Length","metre",
    "mass","kilogram",
    "Mass","kilogram",
    "time","second",
    "Time","second",
    "Temperature","K",
    "tmeperature","K",
    "electric_current","A",
    "Electric_current","A",
    "Amount_of_substance","mole",
    "amount_of_substance","mole",
    "Luminous_intensity","candela",
    "luminous_intensity","candela",
    "area","m^2",
    "density","kg/(m^3)",
    "energy","J",
    "power","W",
    "stress","Pa",
    "capacity","m^3",
    "density","kg/(m^3)",
    "flow","kg/s",
    "force","N",
    "Power","J",
    "power","J",
    "Velocity","m/s",
    "velocity","m/s",
    "speed","m/s",
    "Speed","m/s",
    "viscocityD","Pa*s",
    "viscocityK","(m^2)/s",
    "volume","m^3",
    "Volume","m^3",
    "pressure","Pa",
    "Pressure","Pa",
    "heat","J",
    "Heat","J"
  };

  namespace { // size of tables of units
    const int basic_unit_size = sizeof(UNIT_type::basic_unit_table)/
      sizeof(UNIT_type::basic_units);
 
    const int composite_unit_size = sizeof(UNIT_type::composite_unit_table)/
      sizeof(UNIT_type::composite_units);

    const int cgs_basic_unit_size = sizeof(UNIT_type::cgs_basic_unit_table)/
      sizeof(UNIT_type::basic_units);
 
    const int cgs_composite_unit_size = sizeof(UNIT_type::cgs_composite_unit_table)/sizeof(UNIT_type::composite_units);

    const int reference_unit_size = sizeof(UNIT_type::reference_unit_table)/
      sizeof(UNIT_type::reference_units);

    const int default_unit_size = sizeof(UNIT_type::default_unit_table)/
      sizeof(UNIT_type::default_units);

  } 


  bool UNIT_type::is_reference_unit(string str){
    for(int i=0;i<reference_unit_size;++i){
      if(reference_unit_table[i].name==str)
	return true;
    }
    return false;
  }

  int UNIT_type::where_reference_unit(string str){
    for(int i=0;i<reference_unit_size;++i){
      if(reference_unit_table[i].name==str)
	return i;
    }
    return -1;
  }

  bool UNIT_type::is_composite_unit(string str){
    if(mode==MKS){
      for(int i=0;i<composite_unit_size;++i){
	if(composite_unit_table[i].name==str)
	  return true;
      }
    }
    else if(mode==CGS){
      for(int i=0;i<cgs_composite_unit_size;++i){
	if(cgs_composite_unit_table[i].name==str)
	  return true;
      }
    }
    else
      cout<<"No this mode "<<mode<<endl;
    return false;
  }

  int UNIT_type::where_composite_unit(string str){
    if(mode==MKS){
      for(int i=0;i<composite_unit_size;++i){
	if(composite_unit_table[i].name==str)
	  return i;
      }
    }
    else if(mode==CGS){
      for(int i=0;i<cgs_composite_unit_size;++i){
	if(cgs_composite_unit_table[i].name==str)
	  return i;
      }
    }
    else
      cout<<"No this mode "<<mode<<endl;
    return -1;
  }

  bool UNIT_type::is_basic_unit(string str){
    if(mode==MKS){
      for(int i=0;i<basic_unit_size;++i){
	if(basic_unit_table[i].name==str)
	  return true;
      }
    }
    else if(mode==CGS){
      for(int i=0;i<cgs_basic_unit_size;++i){
	if(cgs_basic_unit_table[i].name==str)
	  return true;
      }
    }
    else
      cout<<"No this mode "<<mode<<endl;
    return false;
  }

  int UNIT_type::where_basic_unit(string str){
    if(mode==MKS){
      for(int i=0;i<basic_unit_size;++i){
	if(basic_unit_table[i].name==str)
	  return i;
      }
    }
    else if(mode==CGS){
      for(int i=0;i<cgs_basic_unit_size;++i){
	if(cgs_basic_unit_table[i].name==str)
	  return i;
      }
    }
    else
      cout<<"No this mode "<<mode<<endl;
    return -1;
  }

  //operation to print the undecomposed unit ----//
  void UNIT_type::printing(ostream &s,exprList exlist)
    {
      if(exlist.size() !=0) {
	exprList::const_iterator i= exlist.begin();
	for (;i!=exlist.end();++i){
	  switch((*i)->op){
	  case OP_NAME:
	    s<<(*i)->name;	
	    {//print "*" between the elements of units
	      if(++i!=exlist.end())
		s<<"*";
	      i--;
	    }
	    break;
	  default:
	    cerr<<"Not the unit operation!"<<endl;
	    break;
	  }
	}
      }
      s<<endl;
    }
  //------build lists of numerator and denominator-----------//
  void UNIT_type::build_lists(exprList &numerator, exprList &denominator, exprP input, bool isnum)
    {
      exprList::const_iterator li,lj;
      switch(input->op){ 
      case OP_NAME:
	if(isnum)
	  numerator.push_back(input);
	else
	  denominator.push_back(input);
	break;
      case OP_TIMES:
	for(li= input->expr_list.begin();li!=input->expr_list.end();++li)
	  build_lists(numerator,denominator,(*li),isnum);
	break;
	case OP_EXOR:
	for(li= input->expr_list.begin();li!=input->expr_list.end();++li){
	  if((*li)->op==OP_INT){
	    //	    cout<<endl<<(*li)->int_val<<endl;
	    lj= input->expr_list.begin();
	    for(int i=1; i!=(*li)->int_val;i++){
	      //cout<<endl<<"i  "<<i<<"  "<<(*lj)->name<<endl;
	      build_lists(numerator,denominator,(*lj),isnum);
	    }
	  }else if((*li)->op==OP_NAME){
	    //cout<<endl<<(*li)->name<<endl;
	    build_lists(numerator,denominator,(*li),isnum);
	  }else
	    cout<<"Wrong exponent!"<<endl;
	}
	break;
      case OP_DIVIDE:
	bool new_isnum=!isnum;
	for(li= input->expr_list.begin();li!=input->expr_list.end();++li){
	  build_lists(numerator,denominator,(*li),isnum);
	  isnum=new_isnum;
	}
	break;
      }
      //cout<<"size of numerator "<<numerator.size()<<endl;
    }

  //count the dimension of numerator or denominator
  void UNIT_type::count_dim(map<string,int> &c_map, const exprList c_list)
    {
      exprList::const_iterator li;
      /*cout<<"content of c_map"<<endl;
      map<string,int>::const_iterator mj=c_map.begin();
      for(;mj!=c_map.end();++mj)
      cout<<mj->first<<"  "<<mj->second<<endl;*/

      for(li= c_list.begin();li!=c_list.end();++li)
	{
	  int dim=1;
	  //cout<<"content of c_list"<<(*li)->name<<endl;
	  if(c_map.find((*li)->name)!=c_map.end())
	    c_map[(*li)->name]+=dim;
	  else
	    c_map[(*li)->name]=1;
	}
      /*map<string,int>::const_iterator mi=c_map.begin();
      for(;mi!=c_map.end();++mi)
      cout<<mi->first<<"  "<<mi->second<<endl;*/
    }

  //remove the duplicated units 
  void UNIT_type::rem_dup(map<string,int> &num_map,map<string,int> &den_map)
    {
      map<string,int>::iterator mi,mj;
      for(mi= num_map.begin();mi!=num_map.end();++mi)
	{
	  for(mj= den_map.begin();mj!=den_map.end();++mj)
	    {
	      if((*mj).first==(*mi).first)
		{
		  int tmp=(*mi).second-(*mj).second;
		  if(tmp>0){
		    num_map[(*mi).first]=tmp;
		    den_map.erase(mj);
		  }
		  else if(tmp<0){
		    den_map[(*mi).first]=-tmp;
		    num_map.erase(mi);
		  }
		  else{
		    num_map.erase(mi);
		    den_map.erase(mj);
		  }
		}
	    }
	}
      /*cout<<"Numerator: "<<endl;
      for(mi= num_map.begin();mi!=num_map.end();++mi)
	cout<<mi->first<<"  "<<mi->second<<endl;
      cout<<"Denominator: "<<endl;
      for(mj= den_map.begin();mj!=den_map.end();++mj)
      cout<<mj->first<<"  "<<mj->second<<endl;*/
    }

  //combine the elements of units//
  map<string,int> UNIT_type::combine_units(map<string,int> num_map,map<string,int> den_map)
    {
      map<string,int>::iterator mi,mj;
      if(num_map.size()!=0){
	for(mi= num_map.begin();mi!=num_map.end();++mi)
	  {
	    for(mj= den_map.begin();mj!=den_map.end();++mj)
	      {
		if((*mj).first==(*mi).first){
		  int tmp=(*mi).second+(*mj).second;
		  if(tmp>0)
		    num_map[(*mi).first]=tmp;
		} else
		  num_map.insert(*mj);
	      }
	  }
      }else{//if there is no element in num_map, insert the elements of den_map
	for(mj= den_map.begin();mj!=den_map.end();++mj){
	  num_map.insert(*mj);
	}
      }
      return num_map;
    }

  //seperate the input units into maps of numerator and denominator//
  void UNIT_type::seperate_unit(map<string,int> &num_map, map<string,int> &den_map,exprP p)    
    {
      exprList numerator,denominator;
      build_lists(numerator,denominator,p);
      //cout<<"Numerator:  ";
      //printing(cout,numerator);
      count_dim(num_map,numerator);
      // cout<<"Denominator:  ";
      //printing(cout,denominator);
      count_dim(den_map,denominator);
      //cout<<endl;
      //cout<<"Remove the duplicated units"<<endl;
      rem_dup(num_map,den_map);
    }

  //change map of numerator or denominator to basic type
  void UNIT_type::change_to_basic_unit(map<string,int>initial_map,map<string,int>&num_map,map<string,int>&den_map,double &conversion_factor)
    {
      map<string,int>::iterator mi;
      exprP exp2;
      exprList numerator, denominator;
      conversion_factor=1;
      //for MKS system
      if(mode==MKS){
	for(mi= initial_map.begin();mi!=initial_map.end();++mi)
	  {
	  if(is_reference_unit((*mi).first)){ 
	    for(int i=0;i!=(*mi).second;i++){//for several same units
	      //cout<<"IS reference unit"<<endl;       
	      int where=where_reference_unit((*mi).first);
	      conversion_factor=conversion_factor*reference_unit_table[where].convert_factor;
	      
	      if(is_composite_unit(reference_unit_table[where].refer_unit)){
		string comp=reference_unit_table[where].refer_unit;
		int where=where_composite_unit(comp);
		conversion_factor=conversion_factor*composite_unit_table[where].convert_factor;
		exp2=expression::create(composite_unit_table[where].derived_unit);
		seperate_unit(num_map,den_map,exp2);
	      }
	      
	      //cout<<"conversion factor:"<<conversion_factor<<endl;
	    }
	  }
	  else if(is_composite_unit((*mi).first)){
	    
	    for(int i=0;i!=(*mi).second;i++){
	      //cout<<"IS composite unit"<<endl;
	      //cout<<(*mi).first<<(*mi).second<<endl;
	      int where=where_composite_unit((*mi).first);
	      conversion_factor=conversion_factor*composite_unit_table[where].convert_factor;
	      exp2=expression::create(composite_unit_table[where].derived_unit);
	      seperate_unit(num_map,den_map,exp2);
	      //cout<<"conversion factor:"<<conversion_factor<<endl;
	    }

	  }
	  else if(is_basic_unit((*mi).first)){
	    //cout<<"IS basic unit"<<endl;
	    if(num_map.find((*mi).first)!=num_map.end())
	      num_map[(*mi).first]=num_map[(*mi).first]+(*mi).second;
	    else
	      num_map[(*mi).first]=(*mi).second;
	    //cout<<(*mi).first<<(*mi).second<<endl;
	    //conversion_factor=conversion_factor*1;
	  }
	  else{
	    unit_error(4,"Not in MKS database! ");
	  }
      }
      }
      // for CGS system
      else if(mode==CGS){
	for(mi= initial_map.begin();mi!=initial_map.end();++mi)
	  {
	  if(is_reference_unit((*mi).first)){ 
	    for(int i=0;i!=(*mi).second;i++){//for several same units
	      //cout<<"IS reference unit"<<endl;       
	      int where=where_reference_unit((*mi).first);
	      conversion_factor=conversion_factor*reference_unit_table[where].convert_factor;
	      
	      if(is_composite_unit(reference_unit_table[where].refer_unit)){
		string comp=reference_unit_table[where].refer_unit;
		int where=where_composite_unit(comp);
		conversion_factor=conversion_factor*cgs_composite_unit_table[where].convert_factor;
		exp2=expression::create(cgs_composite_unit_table[where].derived_unit);
		seperate_unit(num_map,den_map,exp2);
	      }
	      
	      //cout<<"conversion factor:"<<conversion_factor<<endl;
	    }
	  }
	  else if(is_composite_unit((*mi).first)){
	    
	    for(int i=0;i!=(*mi).second;i++){
	      //cout<<"IS composite unit"<<endl;
	      //cout<<(*mi).first<<(*mi).second<<endl;
	      int where=where_composite_unit((*mi).first);
	      conversion_factor=conversion_factor*cgs_composite_unit_table[where].convert_factor;
	      exp2=expression::create(cgs_composite_unit_table[where].derived_unit);
	      seperate_unit(num_map,den_map,exp2);
	      //cout<<"conversion factor:"<<conversion_factor<<endl;
	    }

	  }
	  else if(is_basic_unit((*mi).first)){
	    //cout<<"IS basic unit"<<endl;
	    if(num_map.find((*mi).first)!=num_map.end())
	      num_map[(*mi).first]=num_map[(*mi).first]+(*mi).second;
	    else
	      num_map[(*mi).first]=(*mi).second;
	    //cout<<(*mi).first<<(*mi).second<<endl;
	    conversion_factor=conversion_factor*1;
	  }
	  else{
	    unit_error(5,"Not in CGS database! ");
	    conversion_factor=0;
	  }
	}
	}
    }

  //get the units convert to baisc types//
  void UNIT_type::get_conversion(map<string,int> &num_map, map<string,int> &den_map,double &conversion_factor)
    {
      map<string,int> num_map2,den_map2,num_map3,den_map3;
      double num_conversion_factor,den_conversion_factor;
      int num_size,den_size;
  
      num_size=num_map.size();
      den_size=den_map.size();
      change_to_basic_unit(num_map,num_map2,den_map2,num_conversion_factor);
      //cout<<endl;
      change_to_basic_unit(den_map,num_map3,den_map3,den_conversion_factor);
  
      if(num_size>0&&den_size>0)
	conversion_factor=num_conversion_factor/den_conversion_factor;
      else if(num_size>0&&den_size<=0)
	conversion_factor=num_conversion_factor;
      else if(den_size>0&&num_size<=0)
	conversion_factor=1/den_conversion_factor;
  
      //cout<<endl;
      num_map=combine_units(num_map2,den_map3); 
      den_map=combine_units(den_map2,num_map3); 
    }

  //------unit check---------//
  bool UNIT_type::is_in_db(const std::string &str){
    for(int i=0;i!=basic_unit_size;i++){
      if(str==basic_unit_table[i].name)
	return true;
    }
    for(int i=0;i!=composite_unit_size;i++){
      if(str==composite_unit_table[i].name)
	return true;
    }
    for(int i=0;i!=reference_unit_size;i++){
      if(str==reference_unit_table[i].name)
	return true;
    }
    return false;
  }

  int UNIT_type::in_unit_kind(){
    for(int i=0;i!=default_unit_size;i++){
      if(unit_kind==default_unit_table[i].default_type)
	return i;
    }
    return -1;
  }

  //---if no unit input, set the default unit-------//
  exprP UNIT_type::set_default_unit(exprP &in_exp){
    int tmp=in_unit_kind();
    if(tmp!=-1)
      in_exp=expression::create(default_unit_table[tmp].default_name);
    else{
      unit_error(6,"Unit is not compatible with the kind of unit! ");
    }
    if(unit_kind=="general"){
      unit_error(7,"You must input a unit!");
    }
    return in_exp;
  }

  //------input the expression-------//
  exprP UNIT_type::input(istream &in){
    exprP exp;
    
    cout<<"input the expression"<<endl;
    if(check_unit(in, value)){
      in>>input_unit;
      exp=expression::create(input_unit);
    }
    else
      set_default_unit(exp);
    return exp;
  }

  // check if the unit is single temperature //
  int UNIT_type::is_single_temperature(const exprP input_expr){
    if(input_expr->name=="Fahrenheit"||input_expr->name=="F")
      return 1;
    else if(input_expr->name=="Celsius"||input_expr->name=="C")
      return 2;
    else if(input_expr->name=="Rankine"||input_expr->name=="R")
      return 3;
    return 0;
    }

  //change single temperature to basic type -kelvin- by special calculation
  void UNIT_type::calculate_temperature(exprP &input_expr,double &value){
    seperate_unit(unit_num_map,unit_den_map,input_expr); 
    get_conversion(unit_num_map,unit_den_map,conversion_factor);
    if(conversion_factor>0){
      conversion_factor=1;
      switch(is_single_temperature(input_expr)){
      case 1:
	value=(value+459.67)/1.8;
	break;
      case 2:
	value=value+273.15;
	break;
      case 3:
	value=value/1.8;
	break;
      }
    }else
      cout<<"Not a unit"<<endl;
  }

//change basic type -kelvin- to other temperature by special calculation
  void UNIT_type::reverse_calculate_temperature(exprP &input_expr,double &value){
    seperate_unit(unit_num_map,unit_den_map,input_expr); 
    get_conversion(unit_num_map,unit_den_map,conversion_factor);
    if(conversion_factor>0){
      conversion_factor=1;
      switch(is_single_temperature(input_expr)){
      case 1:
	value=value*1.8-459.67;
	break;
      case 2:
	value=value-273.15;
	break;
      case 3:
	value=value*1.8;
	break;
      }
    }else
      cout<<"Not a unit"<<endl;
  }

  //--- input unit expression and output basic units ---//
  void UNIT_type::output(exprP &in_exp){
      if(is_single_temperature(in_exp)!=0)
	calculate_temperature(in_exp, value);
      else{
	//cout<<endl;
	//cout<<"Loci expression:  ";
	//in_exp->Print(cout);
	
	/*cout<<endl;
	  cout<<endl;
	  cout<<"print the unit:  "<<endl;*/
	seperate_unit(unit_num_map,unit_den_map,in_exp); 
	
	/*cout<<endl;
	  cout<<endl;
	  cout<<"change the unit to basic type:"<<endl;*/
	
	//change the reference type and composite type to the basic type//
	get_conversion(unit_num_map,unit_den_map,conversion_factor);
	if(conversion_factor>0){
	  //cout<<endl;
	  //cout<<"Final unit conversion is: "<<endl;
	  //cout<<"Convert  ";
	  //in_exp->Print(cout);
	  //cout<< "  to basic units: "<<endl;
	  rem_dup(unit_num_map,unit_den_map);
	  /*cout<<"conversion factor is: "<<conversion_factor<<endl;
	    cout<<"Input value is: "<<get_value()<<endl;
	    cout<<"final value is: "<<conversion_factor*value<<endl;*/
	}else
	  cout<<"Not an unit!"<<endl;  
      }
  }

  //-- check if there is unit and get the value input--//
  bool UNIT_type::check_unit(istream &in, double &value){

    parse::kill_white_space(in);
  
    if(in.eof()||in.peek()==EOF){
      cout<<"Nothing input"<<endl;
      return false;
    }

    else if(isdigit(in.peek())){
      if(parse::is_real(in)) 
	value=parse::get_real(in);}
    else if(!isdigit(in.peek())){
      value=1;
      cout<<"No input value, set default to 1"<<endl;
    }
    while(!in.eof()&&isspace(in.peek()))//kill white spaces between the value and unit
      in.get();

    input_value=value;
    if(isalpha(in.peek())||parse::is_token(in,"("))
      return true;
    else
      return false;
  }

  //compare the two units, check if they are comparable
  bool UNIT_type::is_compatible(const std::string unit_str){
    UNIT_type sec_unit;
    exprP sec_exp;
    std::map<std::string,int>::iterator fst_mi;
    std::map<std::string,int> fst_n_map=(*this).unit_num_map,fst_d_map=(*this).unit_den_map;
    sec_unit.mode=(*this).mode;
    sec_unit.unit_kind=(*this).unit_kind;
    sec_unit.value=1;

    sec_exp=expression::create(unit_str);
    sec_unit.output(sec_exp);
    //cout<<sec_unit;

    std::map<std::string,int>::iterator sec_mi;
    std::map<std::string,int> sec_n_map=sec_unit.unit_num_map,sec_d_map=sec_unit.unit_den_map;
    bool num_flag=false,den_flag=false;

    for(fst_mi= fst_n_map.begin();fst_mi!=fst_n_map.end();++fst_mi)
      {
	for(sec_mi= sec_n_map.begin();sec_mi!=sec_n_map.end();++sec_mi)
	  {
	    if(((*sec_mi).first==(*fst_mi).first)&&((*sec_mi).second==(*fst_mi).second))
	      num_flag=true;
	  }
      }
    for(fst_mi= fst_d_map.begin();fst_mi!=fst_d_map.end();++fst_mi)
      {
	for(sec_mi= sec_d_map.begin();sec_mi!=sec_d_map.end();++sec_mi)
	  {
	    if(((*sec_mi).first==(*fst_mi).first)&&((*sec_mi).second==(*fst_mi).second))
	      den_flag=true;
	  }
      }
    if(sec_d_map.size()==0&&fst_d_map.size()==0)
      den_flag=true;
    if((num_flag==true)&&(den_flag==true)&&(sec_d_map.size()==fst_d_map.size())&&(sec_n_map.size()==fst_n_map.size()))
      return true;
    return false;
  }

  //compare the input unit with the default unit, check if they are comparable
  bool UNIT_type::private_is_compatible(){
    exprP exp;
    string str;
    if(mode!=check_available){
      set_default_unit(exp);
      str=exp->name;
      if(is_compatible(str))
	return true;
      else
	return false;
    }
    return true;
  }

  double UNIT_type::get_value_in(const std::string unit_str){
    UNIT_type sec_unit;
    exprP sec_exp;

    sec_unit.mode=(*this).mode;
    sec_unit.unit_kind=(*this).unit_kind;
    sec_unit.value=1;

    sec_exp=expression::create(unit_str);
    if(is_single_temperature(sec_exp)!=0){
      reverse_calculate_temperature(sec_exp,value);
      return value;
    }
    else
      sec_unit.output(sec_exp);
    //cout<<sec_unit;

    return value*(*this).conversion_factor/sec_unit.conversion_factor;
  }

}

