###############################################################################
#
# Copyright 2008-2026, Mississippi State University
#
# This file is part of the Loci Framework.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
###############################################################################

include $(LOCI_BASE)/Loci.conf

SOURCE_ROOT = ../../../..
vpath %.cc $(SOURCE_ROOT)/src/FVMAdapt2/library
INCLUDES = -I../../../contrib/doctest -I$(SOURCE_ROOT)/src/include \
           -I$(SOURCE_ROOT)/src/FVMAdapt2/library
LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
OBJS = test_node_transition.o node_transition.o mesh_transfer.o \
       node_remap.o
NODE_HEADERS = $(SOURCE_ROOT)/src/include/FVMAdapt2/node_transition.h \
               $(SOURCE_ROOT)/src/include/FVMAdapt2/mesh_transfer.h

.PHONY: default clean distclean FRC

default: TestResults

TestResults: test_node_transition.exe FRC
	@LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	  ./test_node_transition.exe > serial.log 2>&1 || { \
	    status=$$?; \
	    cat serial.log; \
	    echo "FVMAdapt2Tests/Core/NodeTransition: FAILED!" > $@; \
	    cat $@; \
	    exit $$status; \
	  }
	@LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	  $(MPI_RUN) -np 2 ./test_node_transition.exe > parallel_2.log 2>&1 || { \
	    status=$$?; cat parallel_2.log; exit $$status; \
	  }
	@LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	  $(MPI_RUN) -np 3 ./test_node_transition.exe > parallel_3.log 2>&1 || { \
	    status=$$?; cat parallel_3.log; exit $$status; \
	  }
	@echo "FVMAdapt2Tests/Core/NodeTransition: PASSED!" > $@
	@cat $@

test_node_transition.exe: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LIBS) $(LDFLAGS)

test_node_transition.o: test_node_transition.cc $(NODE_HEADERS)
	$(CXX) $(COPT) $(DEFINES) $(INCLUDES) $(LOCI_INCLUDES) \
	  -o $@ -c $<

node_transition.o: $(SOURCE_ROOT)/src/FVMAdapt2/library/node_transition.cc \
                   $(NODE_HEADERS)
	$(CXX) $(COPT) $(DEFINES) $(INCLUDES) $(LOCI_INCLUDES) \
	  -o $@ -c $<

mesh_transfer.o: $(SOURCE_ROOT)/src/FVMAdapt2/library/mesh_transfer.cc \
                 $(SOURCE_ROOT)/src/include/FVMAdapt2/mesh_transfer.h
	$(CXX) $(COPT) $(DEFINES) $(INCLUDES) $(LOCI_INCLUDES) \
	  -o $@ -c $<

node_remap.o: $(SOURCE_ROOT)/src/FVMAdapt2/library/node_remap.cc \
                   $(SOURCE_ROOT)/src/FVMAdapt2/library/mesh_state.h \
                   $(NODE_HEADERS)
	$(CXX) $(COPT) $(DEFINES) $(INCLUDES) $(LOCI_INCLUDES) \
	  -o $@ -c $<

clean:
	rm -fr *.o *.exe *.d *.log TestResults debug

distclean: clean

FRC:

ifeq ($(filter $(MAKECMDGOALS),clean distclean),)
-include $(OBJS:.o=.d)
endif
