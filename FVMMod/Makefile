include $(LOCI_BASE)/Loci.conf

INCLUDES = 

#COPT=-O0 #-DBUG1 -DBUG2

# Put objects in the module here
OBJS =\
	metrics.o cellgradient.o facegradient.o limiter.o muscl.o\
	node_values.o bcstuff.o surfaceIntegrations.o petsc.o misc.o




MODULE_NAME = fvm

LOCAL_LIBS = 

JUNK = *~  core ti_files ii_files rii_files

LIB_OBJS=$(OBJS:.o=_lo.o)

default: $(MODULE_NAME)_m.so

$(MODULE_NAME)_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)_m.so $(LIB_FLAGS) $(LIB_OBJS)
FRC : 

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)_m.so $(JUNK)

distclean: 
	rm $(DEPEND_FILES)
	rm -fr $(OBJS) $(TARGET) $(JUNK)


DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
include $(DEPEND_FILES)
