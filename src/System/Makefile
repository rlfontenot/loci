#
# System Makefile
# 

TARGET = libLoci.a #libLoci.so

LOCAL_INCLUDES = -I$(LOCI_BASE)/include


OBJS = Loci_version.o variable.o store_rep.o Map.o constraint.o rule.o \
       fact_db.o execute.o depend_graph.o scheduler.o sched_tools.o \
       sched_mlg.o sched_types.o sched_decompose.o sched_comp.o \
       comp_impl.o comp_reduce.o comp_recurse.o comp_internal.o \
       comp_dag.o comp_cond.o comp_loop.o comp_tools.o distribute.o



JUNK = *~ *.ii core ti_files ii_files rii_files

.SUFFIXES: .o .c .cc .f 

.c.o:    ;  $(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.c
.cc.o:   ;  $(CP) $(CFLAGS) $(TEMPLATE_LIB_OPT) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.cc
.f.o:    ;  $(F77) $(FFLAGS) -c $*.f


default:
	cd .. ; make -f Makefile System


compile: $(TARGET)

libLoci.a: $(OBJS)
	$(AR) $(ARFLAGS) libLoci.a $(OBJS) #../Tools/libTools.a
	$(RANLIB) libLoci.a

libLoci.so: $(OBJS)
	$(LD) $(LDFLAGS) -shared -o libLoci.so $(OBJS) ../Tools/libTools.a

scheduler.o: scheduler.cc
	$(CP) $(CFLAGS) $(TEMPLATE_LIB_OPT) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c scheduler.cc

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)


make.depend:
	touch make.depend 

depend: make.depend
	$(MAKEDEPEND) $(CFLAGS) $(TEMPLATE_LIB_OPT) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) *.cc > make.depend

include make.depend				
