#
# 
include sys.conf
include comp.conf

DEFINES = $(DEBUG) $(SYSTEM_DEFINES) $(MACHINE_SPECIFIC) 

LOCI_BASE = ../

INCLUDES =  

CP  = $(CPP)

COPT     = $(CC_LIB_FLAGS) $(CC_DEBUG) $(CC_OPT1)#$(CC_OPT)
LOPT     = $(COPT)
CFLAGS   = $(COPT) 

F77 = f77
FFLAGS   = -O3 


LIBS     =  $(BASE_LIBS)


default: Tools System hdf5CC

install: default
	LOCI_INSTALL_DIR=$(INSTALL_DIR) bash Install.bash

Tools : FRC
	cd Tools; $(MAKE) compile CC="$(CC)" CP="$(CP)" F77="$(F77)"\
			  DEFINES="$(DEFINES)" INCLUDES="$(INCLUDES)"\
			  CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)"\
			  AR="$(AR)" ARFLAGS="$(ARFLAGS)"\
                          RANLIB="$(RANLIB)" LD="$(LD)"\
			  LDFLAGS="$(LDFLAGS)" LIBS="$(LIBS)" \
			  LOCI_BASE="$(LOCI_BASE)"

hdf5CC : FRC
	cd hdf5CC; $(MAKE) compile CC="$(CC)" CP="$(CP)" F77="$(F77)"\
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
	rm -f hdf5CC/make.depend hdf5CC/make.depend.bak
	rm -f include/*.h~ include/Tools/*.h~
	rm -f Makefile~


clean:
	rm -f $(JUNK)
	cd Tools ; $(MAKE) clean
	cd System ; $(MAKE) clean
	cd hdf5CC ; $(MAKE) clean

Tools_depend: FRC
	cd Tools ; $(MAKE) depend CC="$(CC)" CP="$(CP)" F77="$(F77)" \
			   DEFINES="$(DEFINES)" INCLUDES="$(INCLUDES)"\
			   CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)"\
			   AR="$(AR)" ARFLAGS="$(ARFLAGS)"\
                           RANLIB="$(RANLIB)" LD="$(LD)"\
			   LDFLAGS="$(LDFLAGS)" LIBS="$(LIBS)" \
			   LOCI_BASE="$(LOCI_BASE)"\
                           MAKEDEPEND="$(MAKEDEPEND)"

hdf5CC_depend: FRC
	cd hdf5CC ; $(MAKE) depend CC="$(CC)" CP="$(CP)" F77="$(F77)" \
			   DEFINES="$(DEFINES)" INCLUDES="$(INCLUDES)"\
			   CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)"\
			   AR="$(AR)" ARFLAGS="$(ARFLAGS)"\
                           RANLIB="$(RANLIB)" LD="$(LD)"\
			   LDFLAGS="$(LDFLAGS)" LIBS="$(LIBS)" \
			   LOCI_BASE="$(LOCI_BASE)"\
                           MAKEDEPEND="$(MAKEDEPEND)" 

System_depend: FRC
	cd System ; $(MAKE) depend CC="$(CC)" CP="$(CP)" F77="$(F77)" \
			   DEFINES="$(DEFINES)" INCLUDES="$(INCLUDES)"\
			   CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)"\
			   AR="$(AR)" ARFLAGS="$(ARFLAGS)"\
                           RANLIB="$(RANLIB)" LD="$(LD)"\
			   LDFLAGS="$(LDFLAGS)" LIBS="$(LIBS)" \
			   LOCI_BASE="$(LOCI_BASE)"\
                           MAKEDEPEND="$(MAKEDEPEND)" 

depend: Tools_depend System_depend hdf5CC_depend
