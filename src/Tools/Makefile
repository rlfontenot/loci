#
# Tools Makefile
# 

TARGET = libTools.a #libTools.so



LOCAL_INCLUDES = -I$(LOCI_BASE)/include -I$(LOCI_BASE)/include/hdf5CC

OBJS = fpe.o ftrn_reader.o identity.o eventNotify.o debugger.o \
       parse.o options_list.o expr.o intervalSet.o digraph.o unit_type.o \
	fpe_sun.o

.SUFFIXES: .o .c .cc .f

.c.o:    ;  $(CC) -O $(LIB_FLAGS) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.c
.cc.o:   ;  $(CP) $(LIB_FLAGS) $(CFLAGS) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.cc
.f.o:    ;  $(F77) $(FFLAGS) -c $*.f

default:
	cd .. ; make -f Makefile Tools

compile: $(TARGET)

libTools.a: $(OBJS) 
	$(AR) $(ARFLAGS) libTools.a $(OBJS) 
	$(RANLIB) libTools.a

libTools.so: $(OBJS)
	$(LD) $(LDFLAGS) -shared -o libTools.so $(LIB_FLAGS) $(OBJS) 
	


JUNK =  *~ core ii_files ti_files KCC_files so_locations rii_files

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)


make.depend:
	touch make.depend 

depend: make.depend
	$(MAKEDEPEND) $(CFLAGS) $(TEMPLATE_LIB_OPT) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) *.cc > make.depend

include make.depend				
