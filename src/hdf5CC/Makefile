#
# Tools Makefile
# 

TARGET = libhdf5CC.a

LOCAL_INCLUDES = -I$(LOCI_BASE)/include -I$(LOCI_BASE)/include/hdf5CC

OBJS  =  H5file.o Exception.o H5object.o Group.o Refcount.o Dataspace.o\
	 Dataset.o AbstractDs.o Datatype.o Enumtype.o Atomtype.o Predtype.o\
	 Inttype.o Floatype.o Strtype.o Comptype.o

.SUFFIXES: .o .c .cc .f

.c.o:    ;  cc -O $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.c
.cc.o:   ;  $(CP) $(CFLAGS) $(TEMPLATE_LIB_OPT) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.cc
.f.o:    ;  $(F77) $(FFLAGS) -c $*.f

default:
	cd .. ; make -f Makefile hdf5CC

compile: $(TARGET)

libhdf5CC.a: $(OBJS)
	$(AR) $(ARFLAGS) libhdf5CC.a $(OBJS) 
	$(RANLIB) $(TARGET)

libhdf5CC.so: $(OBJS)
	$(LD) $(LDFLAGS) -shared -o libhdf5CC.so $(OBJS) 


JUNK =  *~ core ii_files ti_files KCC_files so_locations

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)


make.depend:
	touch make.depend 

depend: make.depend
	$(MAKEDEPEND) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) *.cc > make.depend

include make.depend				
