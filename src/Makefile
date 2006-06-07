#
# 
LOCI_BASE = $(shell pwd)

include Loci.conf

#include sys.conf
#include comp.conf

DEFINES = $(DEBUG) $(SYSTEM_DEFINES) $(MACHINE_SPECIFIC) 

default: library lpp FVMtools

library: shared

static: Tools System
shared: Tools/libTools.so System/libLoci.so 

install: default
	INSTALL_DIR=$(INSTALL_DIR) LOCI_INSTALL_DIR=$(LOCI_INSTALL_DIR) bash Install.bash

Tools : FRC
	cd Tools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libTools.a


System : FRC
	cd System; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libLoci.a

Tools/libTools.so: FRC
	cd Tools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libTools.so

System/libLoci.so: FRC
	cd System; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libLoci.so

lpp: library
	cd lpp; $(MAKE) LOCI_BASE="$(LOCI_BASE)" lpp

FVMtools: library
	cd FVMtools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" default
FRC :

FVMMod: library lpp
	cd FVMMod; $(MAKE) LOCI_BASE="$(LOCI_BASE)" default

JUNK =

clean:
	rm -f $(JUNK)
	cd Tools ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" 
	cd System ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" 
	cd lpp ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)"

distclean:
	rm -f $(JUNK)
	cd Tools ; $(MAKE) distclean LOCI_BASE="$(LOCI_BASE)"
	cd System ; $(MAKE) distclean LOCI_BASE="$(LOCI_BASE)"
	cd lpp ; $(MAKE) distclean LOCI_BASE="$(LOCI_BASE)"


