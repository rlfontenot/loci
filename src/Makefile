#
# 
LOCI_BASE = ../
include sys.conf
include comp.conf

DEFINES = $(DEBUG) $(SYSTEM_DEFINES) $(MACHINE_SPECIFIC) 

default: Tools System

shared: Tools/libTools.so System/libLoci.so 

install: default
	LOCI_INSTALL_DIR=$(INSTALL_DIR) CPP='$(CPP)' bash Install.bash

Tools : FRC
	cd Tools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libTools.a


System : FRC
	cd System; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libLoci.a

Tools/libTools.so: FRC
	cd Tools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libTools.so

System/libLoci.so: FRC
	cd System; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libLoci.so


FRC :


JUNK =

allclean: clean
	rm -f Tools/make.depend Tools/make.depend.bak
	rm -f System/make.depend System/make.depend.bak
	rm -f include/*.h~ include/Tools/*.h~
	rm -f Makefile~


clean:
	rm -f $(JUNK)
	cd Tools ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" 
	cd System ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" 

Tools_depend: FRC
	cd Tools ; $(MAKE) depend  LOCI_BASE="$(LOCI_BASE)" 

System_depend: FRC
	cd System ; $(MAKE) depend  LOCI_BASE="$(LOCI_BASE)" 

depend: Tools_depend System_depend
