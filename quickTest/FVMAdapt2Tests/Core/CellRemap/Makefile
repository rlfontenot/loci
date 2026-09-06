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
INCLUDES = -I../../../contrib/doctest -I$(SOURCE_ROOT)/src/include
LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
OBJS = test_cell_remap.o remap_plan.o
REMAP_HEADERS = $(SOURCE_ROOT)/src/include/FVMAdapt2/remap_plan.h \
                $(SOURCE_ROOT)/src/include/FVMAdapt2/mesh_transfer.h

.PHONY: default clean distclean FRC

default: TestResults

TestResults: test_cell_remap.exe FRC
	@LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	  ./test_cell_remap.exe > test_cell_remap.log 2>&1 || { \
	    status=$$?; \
	    cat test_cell_remap.log; \
	    echo "FVMAdapt2Tests/Core/CellRemap: FAILED!" > $@; \
	    cat $@; \
	    exit $$status; \
	  }
	@echo "FVMAdapt2Tests/Core/CellRemap: PASSED!" > $@
	@cat $@

test_cell_remap.exe: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LIBS) $(LDFLAGS)

test_cell_remap.o: test_cell_remap.cc $(REMAP_HEADERS)
	$(CXX) $(COPT) $(DEFINES) $(INCLUDES) $(LOCI_INCLUDES) \
	  -o $@ -c $<

remap_plan.o: $(SOURCE_ROOT)/src/FVMAdapt2/library/remap_plan.cc \
              $(REMAP_HEADERS)
	$(CXX) $(COPT) $(DEFINES) $(INCLUDES) $(LOCI_INCLUDES) \
	  -o $@ -c $<

clean:
	rm -fr *.o *.exe *.d *.log TestResults debug

distclean: clean

FRC:

ifeq ($(filter $(MAKECMDGOALS),clean distclean),)
-include $(OBJS:.o=.d)
endif
