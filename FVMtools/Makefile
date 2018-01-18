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

OTHER_APPS = plot3d2vog make_periodic ugrid2vog fluent2vog cobalt2vog vogmerge vogcut vogcheck extruder vog2surf ccm2vog refmesh marker refine cgns2vog ugrid2cgns  cgns2ugrid cgns2ensight cgns2surf 

# List your object files here
OBJS =  extract.o extract_ensight.o extract_2dgv.o extract_ascii.o extract_tecplot.o extract_vtk.o extract_fv.o extract_cut.o extract_surf.o extract_mean.o extract_cgns.o

LOCAL_LIBS =

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~)  core 

default:
	$(MAKE) $(TARGET) $(OTHER_APPS)

all: $(TARGET) $(OTHER_APPS)


$(TARGET): $(OBJS) 
	$(LD) -o $(TARGET) $(OBJS) $(TEC360LIB) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
FRC : 

libadf/libadf.a:
	(cd libadf; $(MAKE) libadf.a)

ccm2vog.o: ccm2vog.cc
	   $(CP) $(COPT) $(EXCEPTIONS) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -Ilibadf -c ccm2vog.cc

#hash.o: hash.cc
#	   $(CP) $(COPT) $(EXCEPTIONS) $(DEFINES)  -c hash.cc

hash.o:
	$(CC) $(COPTLESS) -c $(DEFINES)  hash.c

getargs.o:
	$(CC) $(COPTLESS) -c $(DEFINES)  getargs.c

cgns2vog.o: cgns2vog.cc
	   $(CP) $(COPT) $(EXCEPTIONS) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -c cgns2vog.cc

#cgns2ugrid.o: cgns2ugrid.cc
#	   $(CP) $(COPT) $(EXCEPTIONS) $(DEFINES) $(INCLUDES) -c cgns2ugrid.c

cgns2ugrid.o:
	$(CC) $(COPTLESS) -c $(DEFINES) cgns2ugrid.c

cgns2ensight.o: cgns2ensight.cc
	   $(CP) $(COPT) $(EXCEPTIONS) $(DEFINES)  $(INCLUDES) -c cgns2ensight.cc

cgns2surf.o: cgns2surf.cc
	   $(CP) $(COPT) $(EXCEPTIONS) $(DEFINES)  $(INCLUDES) -c cgns2surf.cc

ugrid2cgns.o: ugrid2cgns.cc
	   $(CP) $(COPT) $(EXCEPTIONS) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -c ugrid2cgns.cc

#ugrid2cgns.o:
#	$(CC) $(COPTLESS) -c $(DEFINES) ugrid2cgns.c

extract_cgns.o: extract_cgns.cc
	   $(CP) $(COPT) $(EXCEPTIONS) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -c extract_cgns.cc

extract_tecplot.o: extract_tecplot.cc
	   $(CP) $(COPT) $(EXCEPTIONS) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) $(TEC360INCLUDE) -c extract_tecplot.cc

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
cgns2vog: cgns2vog.o hash.o vogtools.o 
	 $(LD) -o cgns2vog cgns2vog.o hash.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cgns2ugrid: cgns2ugrid.o hash.o binaryio.o getargs.o 
	 $(LD) -o cgns2ugrid cgns2ugrid.o hash.o binaryio.o getargs.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cgns2ensight: cgns2ensight.o
	 $(LD) -o cgns2ensight cgns2ensight.o  $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cgns2surf: cgns2surf.o hash.o
	 $(LD) -o cgns2surf cgns2surf.o hash.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)


ugrid2cgns: ugrid2cgns.o hash.o binaryio.o getargs.o
	 $(LD) -o ugrid2cgns ugrid2cgns.o hash.o binaryio.o getargs.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

predicates.o: predicates.c
	$(CC) $(COPTLESS) -c $(DEFINES) predicates.c

make_periodic.o: make_periodic.cc
	$(CPP) $(COPTLESS) -c $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -c make_periodic.cc

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
