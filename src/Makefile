#
# 
LOCI_BASE = ../
include sys.conf
include comp.conf

DEFINES = $(DEBUG) $(SYSTEM_DEFINES) $(MACHINE_SPECIFIC) 

default: Tools System lpp #Tools/libTools.so System/libLoci.so

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

lpp: FRC
	cd lpp; $(MAKE) LOCI_BASE="$(LOCI_BASE)" lpp

FRC :


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


