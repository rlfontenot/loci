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
    /*abbreviation of SI
    "m",Length,1,
    "kg",Mass,1,
    "s",Time,1,
    "K",Temperature,1,
    "A",Electric_current,1,
    "mol",Amount_of_substance,1,
    "cd",Luminous_intensity,1*/
    
  };

  UNIT_type::composite_units UNIT_type::composite_unit_table[]={
    //abbreviation of SI
    "m","metre",1,
    "kg","kilogram",1,
    "s","second",1,
    "K","kelvin",1,
    "A","ampere",1,
    "mol","mol",1,
    "cd","candela",1,

    //metric system
    "centimetre","metre",0.01,//length
    "kilometre","metre",1000,
    "milimetre","metre",0.001,
    "minute","second",60,//time
    "hour","second",3600,
    "day","second",86400,
    "year","second",31536000,
    "pound","kilogram",0.4535924,
    "lb","kilogram",0.4535924,

    "Fahrenheit","kelvin",0.5555556,//temperature interval
    "F","kelvin",0.5555556,
    "Rankine","kelvin",0.5555556,
    "R","kelvin",0.5555556,
    "Celsius","kelvin",1,
    "C","kelvin",1,

    //UK system
    "feet","metre",0.3048,
    "ft","metre",0.3048,
    "mile","metre",1609.344,
    "yard","metre",0.9144,

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
    "atm","pascal",101325,
    "bar","pascal",100000,
    "kPa","pascal",1000,
    "Btu","joule",1055.87
  };

  namespace { // size of tables of units
    const int basic_unit_size = sizeof(UNIT_type::basic_unit_table)/
      sizeof(UNIT_type::basic_units);
 
    const int composite_unit_size = sizeof(UNIT_type::composite_unit_table)/
      sizeof(UNIT_type::composite_units);

    const int reference_unit_size = sizeof(UNIT_type::reference_unit_table)/
      sizeof(UNIT_type::reference_units);
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
    for(int i=0;i<composite_unit_size;++i){
      if(composite_unit_table[i].name==str)
	return true;
    }
    return false;
  }
  int UNIT_type::where_composite_unit(string str){
    for(int i=0;i<composite_unit_size;++i){
      if(composite_unit_table[i].name==str)
	return i;
    }
    return -1;
  }

  bool UNIT_type::is_basic_unit(string str){
    for(int i=0;i<basic_unit_size;++i){
      if(basic_unit_table[i].name==str)
	return true;
    }
    return false;
  }
  int UNIT_type::where_basic_unit(string str){
    for(int i=0;i<basic_unit_size;++i){
      if(basic_unit_table[i].name==str)
	return i;
    }
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
      exprList::const_iterator li;
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
      //cout<<"Denominator:  ";
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
	    conversion_factor=conversion_factor*1;
	  }
	  else{
	    cout<<"Not an unit!"<<endl;
	    abort();
	    conversion_factor=0;
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

  //------input the expression-------//
  exprP UNIT_type::input(istream &in){
    exprP exp;
    cout<<"input the expression"<<endl;
    if(check_unit(in, value))
      exp=expression::create(in);
    else{
      cout<<"No unit input!"<<endl;
      abort();
      //exp=expression::create(in);  
    }
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

  //--- input unit expression and output basic units ---//
  void UNIT_type::output(exprP &in_exp){
      if(is_single_temperature(in_exp)!=0)
	calculate_temperature(in_exp, value);
      else{
	/*cout<<endl;
	  cout<<"Loci expression:  ";
	  exp1->Print(cout);*/
	
	/*cout<<endl;
	  cout<<endl;
	  cout<<"print the unit:  "<<endl;*/
	seperate_unit(unit_num_map,unit_den_map,in_exp); 
	
	/*cout<<endl;
	  cout<<endl;
	  cout<<"change the unit to basic type:"<<endl;*/
	
	//chnage the reference type and composite type to the basic type//
	get_conversion(unit_num_map,unit_den_map,conversion_factor);
	if(conversion_factor>0){
	  cout<<endl;
	  cout<<"Final unit conversion is: "<<endl;
	  cout<<"Convert  ";
	  in_exp->Print(cout);
	  cout<< "  to basic units: "<<endl;
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

    if(isalpha(in.peek()))
      return true;
    else
      return false;
  }

  //get the unit value//
  double UNIT_type::get_value(){
    return value;
  }
}

