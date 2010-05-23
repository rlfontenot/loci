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

OTHER_APPS = cfd++2xdr cobalt2xdr solidMesh2xdr plot3d2xdr pb xdr2cobalt xdr2vog vog2xdr ugrid2vog fluent2vog cobalt2vog vogmerge vogcheck extruder vog2surf

# List your object files here
OBJS =  extract.o extract_ensight.o extract_2dgv.o extract_ascii.o extract_tecplot.o extract_fv.o extract_cut.o extract_surf.o extract_mean.o



LOCAL_LIBS =

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~)  core 

default:
	$(MAKE) $(TARGET) $(OTHER_APPS)

all: $(TARGET) $(OTHER_APPS)


$(TARGET): $(OBJS) 
	$(LD) -o $(TARGET) $(OBJS) $(TEC_LIBS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
FRC : 

extruder: extruder.o
	$(LD) -o extruder extruder.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vogcheck: vogCheck.o
	$(LD) -o vogcheck vogCheck.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

xdr2cobalt: xdr2cobalt.o
	$(LD) -o xdr2cobalt xdr2cobalt.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

xdr2vog: xdr2vog.o vogtools.o
	$(LD) -o xdr2vog xdr2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vog2surf: vog2surf.o vogtools.o
	$(LD) -o vog2surf vog2surf.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

ugrid2vog: ugrid2vog.o vogtools.o
	$(LD) -o ugrid2vog ugrid2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cobalt2vog: cobalt2vog.o vogtools.o
	$(LD) -o cobalt2vog cobalt2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

fluent2vog: fluent2vog.o vogtools.o
	$(LD) -o fluent2vog fluent2vog.o vogtools.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vog2xdr: vog2xdr.o
	$(LD) -o vog2xdr vog2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

vogmerge: vogmerge.o
	$(LD) -o vogmerge vogmerge.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cobalt2xdr: cobalt2xdr.o
	$(LD) -o cobalt2xdr cobalt2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

cfd++2xdr: cfd++2xdr.o
	$(LD) -o cfd++2xdr cfd++2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

plot3d2xdr: plot3d2xdr.o
	$(LD) -o plot3d2xdr plot3d2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

solidMesh2xdr: solidMesh2xdr.o
	$(LD) -o solidMesh2xdr solidMesh2xdr.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

pb: pb.o predicates.o
	$(LD) -o pb pb.o predicates.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

predicates.o:
	$(CC) $(COPT) -c $(DEFINES) predicates.c

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(OTHER_APPS) $(JUNK) *.o

distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES) $(OTHER_APPS) *.o

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
