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

default: library lpp FVMtools FVMMod FVMAdapt FVMOverset

library: shared

static: Tools System
shared: Tools/libTools.so System/libLoci.so METIS SPRNG 

docs:
	cd Tutorial/docs; $(MAKE) LOCI_BASE="$(LOCI_BASE)"

install: default docs
	INSTALL_DIR=$(INSTALL_DIR) LOCI_INSTALL_DIR=$(LOCI_INSTALL_DIR) INSTALL_METIS=$(INSTALL_METIS) bash Install.bash

Tools : FRC
	$(MAKE) -C Tools LOCI_BASE="$(LOCI_BASE)" libTools.a


System : FRC
	$(MAKE) -C System LOCI_BASE="$(LOCI_BASE)" libLoci.a

METIS : FRC
ifeq ($(INSTALL_METIS),1)
	$(MAKE) -C ParMetis-4.0 LOCI_BASE="$(LOCI_BASE)"
endif

SPRNG : FRC
	$(MAKE) -C sprng LOCI_BASE="$(LOCI_BASE)"

Tools/libTools.so: FRC
	$(MAKE) -C Tools LOCI_BASE="$(LOCI_BASE)" libTools.so

System/libLoci.so: FRC
	$(MAKE) -C System LOCI_BASE="$(LOCI_BASE)" libLoci.so

lpp: library
	$(MAKE) -C lpp LOCI_BASE="$(LOCI_BASE)" lpp

FVMtools: library lpp FVMAdapt 
	$(MAKE) -C FVMtools LOCI_BASE="$(LOCI_BASE)" default
FRC :

FVMMod: library lpp
	cd FVMMod; $(MAKE) LOCI_BASE="$(LOCI_BASE)" default

FVMAdapt: library
	cd FVMAdapt; $(MAKE) LOCI_BASE="$(LOCI_BASE)" all

FVMOverset: library lpp
	cd FVMOverset; $(MAKE) LOCI_BASE="$(LOCI_BASE)" all

JUNK =

clean:
	rm -f $(JUNK)
	$(MAKE) -C Tools clean LOCI_BASE="$(LOCI_BASE)" 
	$(MAKE) -C System clean LOCI_BASE="$(LOCI_BASE)" 
	$(MAKE) -C ParMetis-4.0 clean LOCI_BASE="$(LOCI_BASE)" 
	$(MAKE) -C sprng clean LOCI_BASE="$(LOCI_BASE)" 
	$(MAKE) -C FVMtools clean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMMod clean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMAdapt clean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMOverset clean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C lpp clean LOCI_BASE="$(LOCI_BASE)"

distclean:
	rm -f $(JUNK)
	$(MAKE) -C Tools distclean LOCI_BASE="$(LOCI_BASE)" 
	$(MAKE) -C System distclean LOCI_BASE="$(LOCI_BASE)" 
	$(MAKE) -C ParMetis-4.0 distclean LOCI_BASE="$(LOCI_BASE)" 
	$(MAKE) -C sprng distclean LOCI_BASE="$(LOCI_BASE)" 
	$(MAKE) -C FVMtools distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMMod distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMAdapt distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C FVMOverset distclean LOCI_BASE="$(LOCI_BASE)"
	$(MAKE) -C lpp distclean LOCI_BASE="$(LOCI_BASE)"

