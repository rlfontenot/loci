//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#include "lpp.h"
#include <unistd.h>
#include <memory>
using namespace std ;

list<string> include_dirs ;

bool prettyOutput = false ;

namespace {
  const char *revision_name = "$Name: rel-4-0-patches $" ;
  
  std::string version() {
    const char *p = revision_name;
    while(*p!=':' && *p!='\0')
      ++p ;
    if(*p!= '\0')
      ++p ;
    while(*p!=' ' && *p!='\0')
      ++p ;
    if(*p!= '\0')
      ++p ;
    std::string rn ;
    while(*p!='$' &&  *p!=' ' && *p!='\0') 
      rn += *p++ ;

    rn += " lpp compiled at " ;
    rn += __DATE__ ;
    rn += " " ;
    rn += __TIME__ ;
    return rn ;
  }
}

void Usage(int argc, char *argv[]) {
  cerr << "Loci Pre-Processor Usage:" << endl ;
  cerr << argv[0] <<" -I<dir> filename.loci  -o filename1.cc<,filename2.cuda>... <-d cpu<,cuda>...>  <-p>  " << endl ;
  exit(-1) ;
}

void split(string& str, char c, vector<string>& str_list){
  std::size_t found=0, pos=0 ;
  if(str.empty())return;
  while( found != std::string::npos){
    found = str.find(c, pos);
    str_list.push_back(str.substr(pos, found));
    pos = found+1;
  }
  return;
}

                       

        
int main(int argc, char *argv[]) {
  bool file_given = false; 
  string filename ;
  bool out_given = false ;
  vector<string> outfiles ;
  vector<string> device_types ;
  string device_list;
  string out_list;
  bool device_given = false;
  vector<pair<string, std::ostream*> > ofs;
  bool oneOutput = true;

  
  
  for(int i=1;i<argc;++i) {
    if(argv[i][0] == '-') {
      if(argv[i][1] == 'I') {
        string dir = &argv[i][2] ;
        include_dirs.push_back(dir) ;
      }
      if(argv[i][1] == 'p') {
        prettyOutput = true ;
      }

      if(argv[i][1] == 'd') {
        if(i+1>argc ||  device_given)
          Usage(argc,argv) ;
        device_list = argv[i+1] ;
        i++;
        device_given = true;
      }
            
      if(argv[i][1] == 'o') {
        if(i+1>argc || out_given)
          Usage(argc,argv) ;
        out_list = argv[i+1] ;
        i++ ;
        out_given = true ;
      }
      if(argv[i][1] == 'v') {
        cout << "Loci version: " << version() << endl ;
      }
      if(argv[i][1] == 'V') {
        cout << "Loci version: " << version() << endl ;
      }
    } else {
      if(file_given == true) {
        cerr << "multiple filenames given" << endl ;
        Usage(argc,argv) ;
      }
      filename = argv[i] ;
      file_given = true ;
    }
  }
  if(!file_given) {
    cerr << "no filename" << endl ;
    Usage(argc,argv) ;
  }

  if(device_given){
    int nsz = device_list.size() ;
    for(int i=0;i<nsz;++i)
      if(device_list[i] >= 'A' || device_list[i] <= 'Z')
        device_list[i] = std::tolower(device_list[i]) ;
    
    split(device_list, ',', device_types);
  }
 
  if(out_given){
    split(out_list, ',', outfiles);
  }
  
  int num_devs = device_types.size();
  int num_outs = outfiles.size();
 
  if(device_given && num_outs > num_devs){
    cerr << "number of output filenames is greater than number of devices " << endl;
    Usage(argc,argv) ;
  }

  if(num_outs > 1 && num_devs != num_outs){
    cerr << "if multiple output filenams given, number of output files should be equal to number of devices " << endl;
    Usage(argc,argv) ;
  }
  if(num_outs > 1) oneOutput = false;
  
     
  
  parseFile parser ;
 
  try {
    
    //setup ofs 
    if(out_given) {
      if(device_given){
        if(num_outs  == num_devs ){
          for(int i = 0; i < num_outs; i++){
            ofstream* file = new ofstream;
            file->open(outfiles[i].c_str(),ios::out) ;
            if(file->fail()) {
              cerr << "unable to open file " << outfiles[i] << " for writing!" << endl ;
              exit(-1) ;
            }
            ofs.push_back(pair<string, ostream*>(device_types[i],file));
          }
        }else if(num_outs > num_devs){
          cerr << "number of output filenams is greater than number of devices " << endl;
          exit(-1);
        }else if(num_outs == 1){
          ofstream* file = new ofstream;
          file->open(outfiles[0].c_str(),ios::out) ;
          if(file->fail()) {
            cerr << "unable to open file " << outfiles[0] << " for writing!" << endl ;
            exit(-1) ;
          }
          for(int i = 0; i < num_devs; i++){
            ofs.push_back(pair<string, ostream*>(device_types[i],file));
          }
        }else{
          cerr << "if multiple output filenams given, number of output files should be equal to number of devices " << endl;
          exit(-1);
        }
      
      } else { //device not given
        if(num_outs > 1){
          cerr << "number of output filenams is greater than number of devices" << endl;
          exit(-1);
        }else{
          ofstream *file = new ofstream;
          file->open(outfiles[0].c_str(),ios::out) ;
          if(file->fail()) {
            cerr << "unable to open file " << outfiles[0] << " for writing!" << endl ;
            exit(-1) ;
          } 
          ofs.push_back(pair<string, ostream*>("cpu",file));
        }
      }
     
    }else{
      if(device_given){ 
        for(int i = 0; i < num_devs; i++){
          ofs.push_back(pair<string, ostream*>(device_types[i],&cout));
        }
      }else{
        ofs.push_back(pair<string, ostream*>(string("cpu"),&cout));
      }
    }

    //process file
    parser.processFile(filename, ofs, oneOutput) ;
  } catch(parseError pe) {
    if(pe.error_type != "syntax error")
      cerr << pe.error_type << endl ;
    if(out_given){
      for(int i = 0; i < num_outs; i++){
        unlink(outfiles[i].c_str()) ;
      }
    }
    exit(-1) ;
  } catch(...) {
    cerr << "Unknown exception caught!" << endl ;
    if(out_given){
      for(int i = 0; i < num_outs; i++){
        unlink(outfiles[i].c_str()) ;
      }
    }
    exit(-1) ;
  }

  //cleanup
  if(oneOutput){
    if(ofs[0].second != &cout){
      delete ofs[0].second;
    }
  }else{
    for(int i  = 0; i < ofs.size(); i++){
      if(ofs[i].second != &cout){
        delete ofs[i].second;
      }
    }
  }
      
  return 0 ;
}
