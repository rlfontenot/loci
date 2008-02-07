//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <Loci.h>
#include <vector>
#include "defines.h"
using std::cout;
using std::endl;
using std::string;
using std::ifstream;







class init_num_original_nodes : public unit_rule{
  param<int> num_original_nodes;
public:
  init_num_original_nodes(){
    name_store("num_original_nodes", num_original_nodes);
    output("num_original_nodes");
    constraint("UNIVERSE");
    
  }
  //parameter, no loop, 
  virtual void compute(const sequence &seq){
    
    
    *num_original_nodes = 0;
  }
}; 
register_rule<init_num_original_nodes> register_init_num_original_nodes;

class apply_num_original_nodes : public apply_rule<param<int>, Loci::Summation<int> >{
  param<int> num_original_nodes;
  const_store<vect3d> pos;
public:
  apply_num_original_nodes(){
    name_store("num_original_nodes", num_original_nodes);
    name_store("pos", pos);
    input("pos");
    input("num_original_nodes");
    output("num_original_nodes");
   
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity cc){
    *num_original_nodes +=1;
  }
}; 
register_rule<apply_num_original_nodes> register_apply_num_original_nodes;






class get_cellPlan : public pointwise_rule{
  const_param<std::string> planfile_par;
  store<std::vector<char> > cellPlan;
  
  public:
  get_cellPlan(){
    name_store("planfile_par", planfile_par);
    name_store("cellPlan", cellPlan);
    
    input("planfile_par");
    output("cellPlan");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
   
   
    
       hid_t file_id;
       entitySet dom = entitySet(seq);
       file_id = H5Fopen((*planfile_par).c_str(), H5F_ACC_RDONLY,
                         H5P_DEFAULT);
       
       Loci::readContainer(file_id,"cellPlan",cellPlan.Rep(),dom) ;
       H5Fclose(file_id);
       do_loop(seq, this); 
      
  }
  void calculate(Entity cc){
    if(cellPlan[cc].size() == 1 && cellPlan[cc][0] == 'C') cellPlan[cc].resize(0);
  } 
};
register_rule<get_cellPlan> register_get_cellPlan;

class process_plan : public pointwise_rule {
  const_store<std::vector<char> > balancedCellPlan ;
  store<std::vector<char> > tmpPlan;
 public: 
  process_plan(){
    name_store("balancedCellPlan", balancedCellPlan);
    name_store("tmpPlan", tmpPlan);
    input("balancedCellPlan");
    output("tmpPlan");
  }
  virtual void compute(const sequence &seq) {
    do_loop(seq, this);
  }
  
  void calculate(Entity cc){
    if(balancedCellPlan[cc].size() == 0) tmpPlan[cc].push_back('C');
    else tmpPlan[cc] = balancedCellPlan[cc];
  }
    

};
register_rule<process_plan> register_process_plan;
  
class cellplan_output_file : public pointwise_rule {
  const_store<std::vector<char> > tmpPlan ;
  const_param<string> outfile_par ;
  store<bool> cellplan_output ;
  
public:
  cellplan_output_file(){
    name_store("tmpPlan", tmpPlan);
    name_store("outfile_par", outfile_par);
    name_store("cellplan_output", cellplan_output);
    input("(tmpPlan,outfile_par)");
    output("cellplan_output");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    
    hid_t file_id;
    
    file_id = Loci::hdf5CreateFile((*outfile_par).c_str(),H5F_ACC_TRUNC,
                                   H5P_DEFAULT, H5P_DEFAULT) ;
    
    if(Loci::MPI_rank == 0)cout << "writing cellplan into " << *outfile_par << endl;
    Loci::writeContainer(file_id,"cellPlan",tmpPlan.Rep()) ;
    Loci::hdf5CloseFile(file_id) ;
    if(Loci::MPI_rank == 0) cout << "Finish writing cellPlan " << endl;

  
    
  }
 
} ;
register_rule<cellplan_output_file> register_cellplan_output_file;
