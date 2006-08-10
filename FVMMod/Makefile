# This is the Loci base directory
LOCI_BASE ?= /usr/tmp/lush/Loci
# Put objects in the module here
OBJS =	metrics.o cellgradient.o facegradient.o limiter.o muscl.o\
	node_values.o bcstuff.o surfaceIntegrations.o petsc.o misc.o
MODULE_NAME = fvm

#############################################################################
# No changes needed below this line

include $(LOCI_BASE)/Loci.conf

JUNK = *~  core ti_files ii_files rii_files

LIB_OBJS=$(OBJS:.o=_lo.o)

default: $(MODULE_NAME)_m.so

$(MODULE_NAME)_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)_m.so $(LIB_FLAGS) $(LIB_OBJS)
FRC : 

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)_m.so 

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~)  core 
# ".cc" files created from .loci files
LOCI_INTERMEDIATE_FILES = $(subst .loci,.cc, $(wildcard *.loci) )

distclean: 
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)_m.so $(JUNK) $(LOCI_INTERMEDIATE_FILES) $(DEPEND_FILES)


DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
include $(DEPEND_FILES)
