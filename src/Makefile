#
# 

include sys.conf
include comp.conf

DEFINES = $(DEBUG) $(MACHINE_SPECIFIC) 

LOCI_BASE = ../

INCLUDES = 

CP  = $(CPP)
CC  = $(CPP)
COPT     = $(CC_LIB_FLAGS) $(CC_DEBUG) $(CC_OPT1)#$(CC_OPT)
LOPT     = $(COPT)
CFLAGS   = $(COPT) 

F77 = f77
FFLAGS   = -O3 


LIBS     =  $(BASE_LIBS)


default: Tools System

Tools : FRC
	cd Tools; $(MAKE) compile CC="$(CC)" CP="$(CP)" F77="$(F77)"\
			  DEFINES="$(DEFINES)" INCLUDES="$(INCLUDES)"\
			  CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)"\
			  AR="$(AR)" ARFLAGS="$(ARFLAGS)"\
                          RANLIB="$(RANLIB)" LD="$(LD)"\
			  LDFLAGS="$(LDFLAGS)" LIBS="$(LIBS)" \
			  LOCI_BASE="$(LOCI_BASE)"

System : FRC
	cd System; $(MAKE) compile CC="$(CC)" CP="$(CP)" F77="$(F77)" \
			   DEFINES="$(DEFINES)" INCLUDES="$(INCLUDES)"\
			   CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)"\
			   AR="$(AR)" ARFLAGS="$(ARFLAGS)"\
                           RANLIB="$(RANLIB)" LD="$(LD)"\
			   LDFLAGS="$(LDFLAGS)" LIBS="$(LIBS)" \
			   LOCI_BASE="$(LOCI_BASE)"

FRC :


JUNK =

allclean: clean
	rm -f Tools/make.depend Tools/make.depend.bak
	rm -f System/make.depend System/make.depend.bak
	rm -f include/*.h~ include/Tools/*.h~
	rm -f Makefile~


clean:
	rm -f $(JUNK)
	cd Tools ; $(MAKE) clean
	cd System ; $(MAKE) clean

Tools_depend: FRC
	cd Tools ; $(MAKE) depend INCLUDES="$(INCLUDES)"\
			  DEFINES="$(DEFINES)" \
                          CC_INCLUDE_PATH="$(CC_INCLUDE_PATH)"\
			  LOCI_BASE="$(LOCI_BASE)"

System_depend: FRC
	cd System ; $(MAKE) depend INCLUDES="$(INCLUDES)"\
			  DEFINES="$(DEFINES)" \
                          CC_INCLUDE_PATH="$(CC_INCLUDE_PATH)"\
			  LOCI_BASE="$(LOCI_BASE)"

depend: Tools_depend System_depend 
