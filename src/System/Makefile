#
# System Makefile
# 

include $(LOCI_BASE)/Loci.conf

TARGET = libLoci.a #libLoci.so

OBJS = Loci_version.o Loci_Datatypes.o variable.o store_rep.o Map.o \
       constraint.o rule.o fact_db.o execute.o depend_graph.o \
       scheduler.o sched_tools.o sched_mlg.o sched_types.o \
       sched_decompose.o sched_comp.o comp_impl.o comp_reduce.o \
       comp_recurse.o comp_internal.o comp_dag.o comp_cond.o \
       comp_loop.o comp_tools.o distribute.o DMap.o DMultiMap.o\
       param_rule.o sched_db.o 

COPT = $(CC_DEBUG) $(CC_LIB_FLAGS) $(CC_OPT1)

default:
	cd .. ; make -f Makefile System

compile: $(TARGET)

libLoci.a: $(OBJS)
	$(AR) $(ARFLAGS) libLoci.a $(OBJS) 
	$(RANLIB) libLoci.a
LIB_OBJS=$(OBJS:.o=_lo.o)
libLoci.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) libLoci.so $(LIB_FLAGS) $(LIB_OBJS) 

JUNK =  *~ core ii_files ti_files KCC_files so_locations rii_files

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK)


make.depend:
	touch make.depend 

depend: make.depend
	$(MAKEDEPEND) $(CFLAGS)  $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) $(LOCAL_INCLUDES) *.cc > make.depend

include make.depend				
