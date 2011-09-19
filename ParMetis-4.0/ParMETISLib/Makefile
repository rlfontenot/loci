
include $(LOCI_BASE)/Loci.conf

INCLUDES=-I. -I.. -I../GKLib -I../METISLib  -DNDEBUG -DNDEBUG2

COPT = $(CC_LIB_FLAGS) $(COPTLESS)

OBJS = \
	akwayfm.o ametis.o balancemylink.o comm.o csrmatch.o ctrl.o \
	debug.o diffutil.o frename.o gkmetis.o gkmpi.o graph.o \
	initbalance.o initmsection.o initpart.o kmetis.o kwayrefine.o \
	match.o mdiffusion.o mesh.o mmetis.o move.o msetup.o node_refine.o \
	ometis.o pspases.o redomylink.o remap.o renumber.o rmetis.o selectq.o \
	serial.o stat.o timer.o util.o wave.o weird.o wspace.o xyzpart.o

libparmetis.a: $(OBJS) 
	$(AR) $(ARFLAGS) libparmetis.a $(OBJS) 
	$(RANLIB) libparmetis.a

LIB_OBJS=$(OBJS:.o=_lo.o)
libparmetis.so: $(LIB_OBJS)
	$(DYNAMIC_LD) $(DYNAMIC_LD_FLAGS) libparmetis.$(LIB_SUFFIX) $(LIB_FLAGS) $(LIB_OBJS) 

JUNK =  *~ core ii_files ti_files KCC_files so_locations rii_files

clean:
	rm -fr $(LIB_OBJS) $(OBJS) $(TARGET) $(JUNK) libparmetis.$(LIB_SUFFIX)


distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES) libparmetis.$(LIB_SUFFIX)

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
