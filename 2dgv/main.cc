#include "Xdev.h"
#include <stdlib.h>
using std::cout ;
using std::cerr ;
using std::endl ;

char *progname ;

void showhelp()
{
  cout << "Usage:" << endl
       << progname << " <options> filename" << endl
       << "Options are:" << endl
       << "-c #" << endl
       << "     This option sets the jump between contour lines." << endl
       << endl
       << "-s" << endl
       << "     This option turns off the default shaded output to Postscript devices."
       << endl << endl
       << "-S #" << endl 
       << "     This option sets the darkness level of output shading."
       << endl << endl
       << "-f #" << endl
       << "     This option sets the label font size for Postscript devices."
       << endl << endl
       << "-v" << endl
       <<"      Read in function values at cell centers (default is vertici centers.)"
       << endl << endl
       << "-d" << endl
       << "     Only read in and display the grid." << endl << endl
       << "-r" << endl 
       << "     Turn Ruler on." << endl
       << "-o hpgl\n-o postscript\n-o printer" << endl
       << "     Send output to either gridout.hpgl, gridout.ps, or postscript printer."
       << endl << endl ;
  exit(-1) ; 
}

int main(int argc, char *argv[]) {
  double font_size ;
  
  enum {XWINDOWS,HPGLFILE,PSFILE,PSPRINTER} device = XWINDOWS ;
  device = XWINDOWS ;
  progname = argv[0] ;
  int j= 0 ;
  for(int i=1; i < argc; i++) {
    if ((strcmp("-f", argv[i]) == 0) && (i+1 <= argc)) {
      i++; j=j+2;
      font_size = atof(argv[i]) ;
      if(font_size < 3 || font_size > 30) {
        cerr << "font size must be between 3 and 30" << endl ;
        font_size = 6 ;
        break ;
      }
    } else if ((strcmp("-o",argv[i]) == 0) && (i+1 <= argc)) {
      i++;
      j=j+2;
      if(!strcmp(argv[i],"hpgl"))
        device = HPGLFILE ;
      else if(!strcmp(argv[i],"postscript"))
        device = PSFILE ;
      else if(!strcmp(argv[i],"printer"))
        device = PSPRINTER ;
      else
        showhelp() ;
    } else if (((!strcmp("-h",argv[i])) || (!strcmp("-help",argv[i])))
               && (i+1 <= argc)) showhelp();
  }
  if(j > argc-2 ) {
    if(device != XWINDOWS)
      showhelp();
  } else {
    grids.push_back(grid()) ;
    grids.front().input(argv[argc-1]) ;
    //    strcpy(gridfilename,argv[argc-1]) ;
    //    readgrid(gridfilename);
  }
  if(device == XWINDOWS)
    xinteract(argc, argv) ;
  //  else
  //    Plotmesh(world,device) ;
  exit(0) ;
}

