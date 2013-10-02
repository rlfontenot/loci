###############################################################################
#
# Copyright 2008, Mississippi State University
#
# This file is part of the Loci Framework.
#
# The Loci Framework is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The Loci Framework is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
#
###############################################################################
#
LOCI_BASE = $(shell pwd)
LOCI_COMPILE = true
include Loci.conf

#include sys.conf
#include comp.conf

DEFINES = $(DEBUG) $(SYSTEM_DEFINES) $(MACHINE_SPECIFIC) 

default: library lpp FVMtools FVMMod FVMAdapt 

library: shared

static: Tools System
shared: Tools/libTools.so System/libLoci.so METIS SPRNG

docs:
	cd Tutorial/docs; $(MAKE) LOCI_BASE="$(LOCI_BASE)"

install: default docs
	INSTALL_DIR=$(INSTALL_DIR) LOCI_INSTALL_DIR=$(LOCI_INSTALL_DIR) bash Install.bash

Tools : FRC
	cd Tools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libTools.a


System : FRC
	cd System; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libLoci.a

METIS : FRC
	cd ParMetis-4.0; $(MAKE) LOCI_BASE="$(LOCI_BASE)" 

SPRNG : FRC
	cd sprng; $(MAKE) LOCI_BASE="$(LOCI_BASE)"

Tools/libTools.so: FRC
	cd Tools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libTools.so

System/libLoci.so: FRC
	cd System; $(MAKE) LOCI_BASE="$(LOCI_BASE)" libLoci.so

lpp: library
	cd lpp; $(MAKE) LOCI_BASE="$(LOCI_BASE)" lpp

FVMtools: library lpp
	cd FVMtools; $(MAKE) LOCI_BASE="$(LOCI_BASE)" default
FRC :

FVMMod: library lpp
	cd FVMMod; $(MAKE) LOCI_BASE="$(LOCI_BASE)" default

FVMAdapt: library
	cd FVMAdapt; $(MAKE) LOCI_BASE="$(LOCI_BASE)" all

JUNK =

clean:
	rm -f $(JUNK)
	cd Tools ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" 
	cd System ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" 
	cd ParMetis-4.0 ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" 
	cd sprng ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" 
	cd FVMtools; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)"
	cd FVMMod; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)"
	cd FVMAdapt; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)"
	cd lpp ; $(MAKE) clean LOCI_BASE="$(LOCI_BASE)"

distclean:
	rm -f $(JUNK)
	cd Tools ; $(MAKE) distclean LOCI_BASE="$(LOCI_BASE)"
	cd System ; $(MAKE) distclean LOCI_BASE="$(LOCI_BASE)"
	cd ParMetis-4.0 ; $(MAKE) distclean LOCI_BASE="$(LOCI_BASE)" 
	cd sprng ; $(MAKE) distclean LOCI_BASE="$(LOCI_BASE)" 
	cd FVMtools; $(MAKE) distclean LOCI_BASE="$(LOCI_BASE)"
	cd FVMMod; $(MAKE) distclean LOCI_BASE="$(LOCI_BASE)"
	cd FVMAdapt; $(MAKE) distclean LOCI_BASE="$(LOCI_BASE)"
	cd lpp ; $(MAKE) distclean LOCI_BASE="$(LOCI_BASE)"

