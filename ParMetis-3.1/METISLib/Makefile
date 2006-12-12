
include $(LOCI_BASE)/Loci.conf

INCLUDES=-I.

COPT = $(CC_LIB_FLAGS) $(COPTLESS)

OBJS = coarsen.o fm.o initpart.o match.o ccgraph.o \
       pmetis.o pqueue.o refine.o util.o timing.o debug.o \
       bucketsort.o graph.o stat.o kmetis.o kwayrefine.o \
       kwayfm.o balance.o ometis.o srefine.o sfm.o separator.o \
       mincover.o mmd.o mesh.o meshpart.o frename.o fortran.o \
       myqsort.o compress.o parmetis.o estmem.o \
       mpmetis.o mcoarsen.o mmatch.o minitpart.o mbalance.o \
       mutil.o mkmetis.o mkwayrefine.o mkwayfmh.o \
       mrefine2.o minitpart2.o mbalance2.o mfm2.o \
       kvmetis.o kwayvolrefine.o kwayvolfm.o subdomains.o \
       mfm.o memory.o mrefine.o checkgraph.o

libmetis.a: $(OBJS) 
	$(AR) $(ARFLAGS) libmetis.a $(OBJS) 
	$(RANLIB) libmetis.a

LIB_OBJS=$(OBJS:.o=_lo.o)
libmetis.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) libmetis.so $(LIB_FLAGS) $(LIB_OBJS) 

JUNK =  *~ core ii_files ti_files KCC_files so_locations rii_files

clean:
	rm -fr $(LIB_OBJS) $(OBJS) $(TARGET) $(JUNK)


distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES)
	

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
include $(DEPEND_FILES)
