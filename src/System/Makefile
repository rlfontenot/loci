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
# System Makefile
# 


LOCI_COMPILE = true

include $(LOCI_BASE)/Loci.conf

TARGET = libLoci.a #libLoci.so
INCLUDES = -I$(LOCI_BASE)/sprng/include

OBJS = Loci_version.o  Initialize.o Loci_Datatypes.o variable.o \
       store_rep.o Map.o constraint.o DMap.o DMultiMap.o \
       rule.o fact_db.o execute.o depend_graph2.o \
       scheduler.o sched_tools.o sched_mlg.o sched_types.o \
       sched_decompose.o sched_comp.o comp_impl.o comp_reduce.o \
       comp_recurse.o comp_internal.o comp_dag.o comp_cond.o \
       comp_loop.o comp_chomp.o comp_tools.o distribute.o categories.o \
       comp_dynamic.o comp_constraint.o \
       comp_map.o param_rule.o sched_db.o  hdf5_readwrite.o mod_db.o \
       visit_assemble.o visit_deco.o visit_tools.o visit_viz.o\
       visit_chomp.o dynamic_exec.o \
       dist_tools.o distribute_io.o distribute_container.o dist_internal.o \
       FVMGridReader.o FVMGridWriter.o FVMStuff.o kd_tree.o pnn.o orb.o \
       keyspace.o interpolate_tools.o thread.o

COPT = $(CC_LIB_FLAGS) $(COPTLESS) -DLOCI_RPATH=\"$(LOCI_RPATH)\"
#COPT = $(CC_LIB_FLAGS) -g -O0 -DLOCI_RPATH=\"$(LOCI_RPATH)\"

default:
	cd .. ; make -f Makefile System

compile: $(TARGET)

libLoci.a: $(OBJS) 
	$(AR) $(ARFLAGS) libLoci.a $(OBJS) 
	$(RANLIB) libLoci.a

LIB_OBJS=$(OBJS:.o=_lo.o)
libLoci.so: $(LIB_OBJS) 
	$(DYNAMIC_LD) $(DYNAMIC_LD_FLAGS) libLoci.$(LIB_SUFFIX) $(LIB_FLAGS) $(LIB_OBJS) 

Loci_version.o: FORCE
Loci_version_lo.o: FORCE

.PHONY: FORCE
FORCE:

JUNK =  *~ core ii_files ti_files KCC_files so_locations rii_files

.PHONY: clean
clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) libLoci.$(LIB_SUFFIX)

.PHONY: distclean
distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES) libLoci.$(LIB_SUFFIX)


DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
