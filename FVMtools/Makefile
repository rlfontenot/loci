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
TARGET = extract

OTHER_APPS = plot3d2vog make_periodic ugrid2vog fluent2vog cobalt2vog vogmerge vogcut vogcheck extruder vog2surf ccm2vog refmesh marker refine #extract_chem 

# List your object files here
OBJS =  extract.o extract_ensight.o extract_2dgv.o extract_ascii.o extract_tecplot.o extract_vtk.o extract_fv.o extract_cut.o extract_surf.o extract_mean.o 

LOCAL_LIBS =

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~)  core 

default:
	$(MAKE) $(TARGET) $(OTHER_APPS)

all: $(TARGET) $(OTHER_APPS)


$(TARGET): $(OBJS) 
	$(LD) -o $(TARGET) $(OBJS) $(TEC_LIBS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
FRC : 

libadf/libadf.a:
	(cd libadf; $(MAKE) libadf.a)

ccm2vog.o: ccm2vog.cc
	   $(CP) $(COPT) $(EXCEPTIONS) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -Ilibadf -c ccm2vog.cc

extruder: extruder.o
	$(LD) -o extruder extruder.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vogcheck: vogCheck.o
	$(LD) -o vogcheck vogCheck.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)


vog2surf: vog2surf.o vogtools.o
	$(LD) -o vog2surf vog2surf.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

extract_chem: extract_chem.o 
	$(LD) -o extract_chem extract_chem.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

ugrid2vog: ugrid2vog.o vogtools.o
	$(LD) -o ugrid2vog ugrid2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

ugridmirror: ugridmirror.o vogtools.o
	$(LD) -o ugridmirror ugridmirror.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cobalt2vog: cobalt2vog.o vogtools.o
	$(LD) -o cobalt2vog cobalt2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

fluent2vog: fluent2vog.o vogtools.o
	$(LD) -o fluent2vog fluent2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

ccm2vog: ccm2vog.o vogtools.o libadf/libadf.a
	 $(LD) -o ccm2vog ccm2vog.o vogtools.o libadf/libadf.a $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vogmerge: vogmerge.o
	$(LD) -o vogmerge vogmerge.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vogcut: vogcut.o vogtools.o
	$(LD) -o vogcut vogcut.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)


plot3d2vog: plot3d2vog.o vogtools.o
	$(LD) -o plot3d2vog plot3d2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

make_periodic: make_periodic.o predicates.o vogtools.o
	$(LD) -o make_periodic make_periodic.o predicates.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

predicates.o:
	$(CC) $(COPTLESS) -c $(DEFINES) predicates.c

refmesh: refmesh.o
	$(LD) -o  refmesh  refmesh.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS) -lfvmadaptfunc
marker: marker.o
	$(LD) -o  marker  marker.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS) -lfvmadaptfunc
refine: refine.o
	$(LD) -o  refine  refine.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS) -lfvmadaptfunc


clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(OTHER_APPS) $(JUNK) *.o

distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES) $(OTHER_APPS) *.o

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
