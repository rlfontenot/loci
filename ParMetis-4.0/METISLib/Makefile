
include $(LOCI_BASE)/Loci.conf

INCLUDES=-I. -I.. -I../GKLib -DNDEBUG -DNDEBUG2

COPT = $(CC_LIB_FLAGS) $(COPTLESS)

OBJS = \
	auxapi.o balance.o bucketsort.o checkgraph.o coarsen.o compress.o \
	contig.o debug.o fm.o fortran.o frename.o gklib.o graph.o initpart.o \
	kmetis.o kwayfm.o kwayrefine.o mcutil.o mesh.o meshpart.o minconn.o \
	mincover.o mmd.o ometis.o options.o parmetis.o pmetis.o refine.o\
	separator.o sfm.o srefine.o stat.o timing.o util.o wspace.o

libmetis.a: $(OBJS) 
	$(AR) $(ARFLAGS) libmetis.a $(OBJS) 
	$(RANLIB) libmetis.a

LIB_OBJS=$(OBJS:.o=_lo.o)
libmetis.so: $(LIB_OBJS)
	$(DYNAMIC_LD) $(DYNAMIC_LD_FLAGS) libmetis.$(LIB_SUFFIX) $(LIB_FLAGS) $(LIB_OBJS) -L../GKLib -lgk

JUNK =  *~ core ii_files ti_files KCC_files so_locations rii_files

clean:
	rm -fr $(LIB_OBJS) $(OBJS) $(TARGET) $(JUNK) libmetis.$(LIB_SUFFIX)


distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES) libmetis.$(LIB_SUFFIX)

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
