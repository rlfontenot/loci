#
# Tools Makefile
# 

include $(LOCI_BASE)/Loci.conf

INCLUDES = -I$(LOCI_BASE)/include/hdf5CC

COPT = $(CC_DEBUG) $(CC_LIB_FLAGS) $(CC_OPT)

TARGET = libhdf5CC.a #libhdf5CC.so


OBJS  =  H5file.o Exception.o H5object.o Group.o Refcount.o Dataspace.o\
	 Dataset.o AbstractDs.o Datatype.o Enumtype.o Atomtype.o Predtype.o\
	 Inttype.o Floatype.o Strtype.o Comptype.o

default:
	cd .. ; make -f Makefile hdf5CC

compile: $(TARGET)

libhdf5CC.a: $(OBJS)
	$(AR) $(ARFLAGS) libhdf5CC.a $(OBJS) 
	$(RANLIB) libhdf5CC.a

LIB_OBJS=$(OBJS:.o=.lo)
libhdf5CC.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) libhdf5CC.so $(LIB_FLAGS) $(LIB_OBJS) 



JUNK =  *~ core ii_files ti_files KCC_files so_locations rii_files

clean:
	rm -fr $(LIB_OBJS) $(OBJS) $(TARGET) $(JUNK)


make.depend:
	touch make.depend 

depend: make.depend
	$(MAKEDEPEND) $(CFLAGS)  $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) $(LOCAL_INCLUDES) *.cc > make.depend

include make.depend				
