#
# 
LOCI_BASE=../
include sys.conf
include comp.conf

DEFINES = $(DEBUG) $(SYSTEM_DEFINES) $(MACHINE_SPECIFIC) 

default: Tools System hdf5CC

shared: #Tools/libTools.so System/libLoci.so hdf5CC/libhdf5CC.so

install: default
	LOCI_INSTALL_DIR=$(INSTALL_DIR) bash Install.bash

Tools : FRC
	cd Tools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libTools.a

hdf5CC : FRC
	cd hdf5CC; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libhdf5CC.a

System : FRC
	cd System; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libLoci.a

Tools/libTools.so: FRC
	cd Tools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libTools.so

System/libLoci.so: FRC
	cd System; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libLoci.so
	
hdf5CC/libhdf5CC.so : FRC
	cd hdf5CC; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libhdf5CC.so


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
