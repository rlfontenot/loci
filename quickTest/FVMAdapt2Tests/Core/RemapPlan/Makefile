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

INCLUDES = -I../../../contrib/doctest -I$(LOCI_BASE)/include \
           -I../../../../src/FVMAdapt2/library
LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
SOURCES = test_face_transition.cc
OBJS = $(SOURCES:.cc=.o)

.PHONY: default clean distclean FRC

default: TestResults

TestResults: test_face_transition.exe FRC
	@LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	  ./test_face_transition.exe > test_face_transition.log 2>&1 || { \
	    status=$$?; \
	    cat test_face_transition.log; \
	    echo "FVMAdapt2Tests/Core/RemapPlan: FAILED!" > $@; \
	    cat $@; \
	    exit $$status; \
	  }
	@echo "FVMAdapt2Tests/Core/RemapPlan: PASSED!" > $@
	@cat $@

test_face_transition.exe: $(OBJS) \
                          $(LOCI_BASE)/lib/libfvmadapt2func.$(LIB_SUFFIX)
	$(LD) -o $@ $(OBJS) -L$(LOCI_BASE)/lib \
	  -lfvmadapt2func $(LIBS) $(LDFLAGS)

clean:
	rm -fr *.o *.exe *.d *.log TestResults debug

distclean: clean

FRC:

ifeq ($(filter $(MAKECMDGOALS),clean distclean),)
-include $(SOURCES:.cc=.d)
endif
