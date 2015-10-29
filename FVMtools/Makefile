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
include $(LOCI_BASE)/Loci.conf

COPT = $(COPTLESS)
#COPT = -g -O0 #-DDEBUG -DBOUNDS_CHECK

# List the name of your compiled program here
#TARGET = refine

OTHER_APPS = refine refmesh marker 

# List your object files here
#OBJS =  extract.o extract_ensight.o extract_2dgv.o extract_ascii.o extract_tecplot.o extract_vtk.o extract_fv.o extract_cut.o extract_surf.o extract_mean.o 

LOCAL_LIBS =

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~)  core 

default:
	$(MAKE)  $(OTHER_APPS)

all:  $(OTHER_APPS)



FRC : 

libadf/libadf.a:
	(cd libadf; $(MAKE) libadf.a)


refmesh: refmesh.o
	$(LD) -o  refmesh  refmesh.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS) -lfvmadaptfunc
marker: marker.o
	$(LD) -o  marker  marker.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS) -lfvmadaptfunc
refine: refine.o
	$(LD) -o  refine  refine.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS) -lfvmadaptfunc


clean:
	rm -fr $(OBJS) $(LIB_OBJS)  $(OTHER_APPS) $(JUNK) *.o

distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES) $(OTHER_APPS) *.o

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
