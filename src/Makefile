#
# 
LOCI_BASE=../
include sys.conf
include comp.conf

DEFINES = $(DEBUG) $(SYSTEM_DEFINES) $(MACHINE_SPECIFIC) 

default: Tools System hdf5CC

install: default
	LOCI_INSTALL_DIR=$(INSTALL_DIR) bash Install.bash

Tools : FRC
	cd Tools; $(MAKE) compile LOCI_BASE="$(LOCI_BASE)"

hdf5CC : FRC
	cd hdf5CC; $(MAKE) compile LOCI_BASE="$(LOCI_BASE)"

System : FRC
	cd System; $(MAKE) compile LOCI_BASE="$(LOCI_BASE)"

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
	cd Tools ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" 
	cd System ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" 
	cd hdf5CC ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)"

Tools_depend: FRC
	cd Tools ; $(MAKE) depend  LOCI_BASE="$(LOCI_BASE)" 

hdf5CC_depend: FRC
	cd hdf5CC ; $(MAKE) depend  LOCI_BASE="$(LOCI_BASE)" 

System_depend: FRC
	cd System ; $(MAKE) depend  LOCI_BASE="$(LOCI_BASE)" 

depend: Tools_depend System_depend hdf5CC_depend
