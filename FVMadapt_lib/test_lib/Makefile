##############################################################################
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
##############################################################################
include $(LOCI_BASE)/Loci.conf

OBJS = refine.o 
MARKEROBJS = marker.o 
REFMESHOBJS = refmesh.o

LIBXML2_LIBS ?= -lxml2
LIBXML2_DEFINE ?= -DUSE_LIBXML2
LIBXML2_INCLUDE ?= -I/usr/include/libxml2
LOCAL_LIBS      = $(LIBXML2_LIBS) -L$(LOCI_BASE)/lib -lfvmadaptfunc
FVMADAPT_INCLUDE ?= -I$(LOCI_BASE)/FVMadapt_lib/include
#LDFLAGS = -Wl -rpath /home/qxue/FVMadapt/

JUNK =		*~ \
		core \
		ti_files \
		ii_files \
		rii_files

INCLUDES=	$(FVMADAPT_INCLUDE) $(LIBXML2_INCLUDE) $(LIBXML2_DEFINE)


default: refine marker refmesh

refine: $(OBJS) $(MODULES)
	$(LD) -o  refine  $(OBJS) $(MODULES) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

marker: $(MARKEROBJS) $(MODULES)
	$(LD) -o  marker  $(MARKEROBJS) $(MODULES) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

refmesh: $(REFMESHOBJS) $(MODULES)
	$(LD) -o  refmesh  $(REFMESHOBJS) $(MODULES) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

FRC :

clean:
	rm -fr $(OBJS) $(MARKEROBJS) $(REFMESHOBJS)  $(TARGET) $(JUNK)

distclean: clean
	rm -fr $(DEPEND_FILES) refine marker refmesh

DEPEND_FILES=$(subst .o,.d,$(OBJS) $(MARKEROBJS) $(REFMRSHOBJS))
#include automatically generated dependencies
-include $(DEPEND_FILES)
