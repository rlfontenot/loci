
include $(LOCI_BASE)/Loci.conf

INCLUDES=-I.

COPT = $(CC_LIB_FLAGS) $(COPTLESS)

OBJS = comm.o util.o debug.o setup.o grsetup.o timer.o \
       node_refine.o initmsection.o order.o \
       xyzpart.o pspases.o frename.o \
       iintsort.o iidxsort.o ikeysort.o ikeyvalsort.o \
       kmetis.o gkmetis.o ometis.o \
       initpart.o match.o coarsen.o \
       kwayfm.o kwayrefine.o kwaybalance.o \
       remap.o stat.o fpqueue.o \
       ametis.o rmetis.o lmatch.o initbalance.o \
       mdiffusion.o diffutil.o wave.o \
       csrmatch.o redomylink.o balancemylink.o \
       selectq.o akwayfm.o serial.o move.o \
       mmetis.o mesh.o memory.o weird.o backcompat.o

libparmetis.a: $(OBJS) 
	$(AR) $(ARFLAGS) libparmetis.a $(OBJS) 
	$(RANLIB) libparmetis.a

LIB_OBJS=$(OBJS:.o=_lo.o)
libparmetis.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) libparmetis.so $(LIB_FLAGS) $(LIB_OBJS) 

JUNK =  *~ core ii_files ti_files KCC_files so_locations rii_files

clean:
	rm -fr $(LIB_OBJS) $(OBJS) $(TARGET) $(JUNK)


distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES)
	

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
