#
# Tools Makefile
# 

include $(LOCI_BASE)/Loci.conf

TARGET = libTools.a libTools.so

COPT = $(CC_DEBUG) $(CC_LIB_FLAGS) $(CC_OPT)


INCLUDES = 

OBJS = fpe.o ftrn_reader.o identity.o eventNotify.o debugger.o \
       parse.o options_list.o expr.o intervalSet.o digraph.o unit_type.o \
	fpe_sun.o

default:
	cd .. ; make -f Makefile Tools

compile: $(TARGET)

libTools.a: $(OBJS) 
	$(AR) $(ARFLAGS) libTools.a $(OBJS) 
	$(RANLIB) libTools.a

LIB_OBJS=$(OBJS:.o=_lo.o)
libTools.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) libTools.so $(LIB_FLAGS) $(LIB_OBJS) 


JUNK =  *~ core ii_files ti_files KCC_files so_locations rii_files

clean:
	rm -fr $(LIB_OBJS) $(OBJS) $(TARGET) $(JUNK)


make.depend:
	touch make.depend 

depend: make.depend
	$(MAKEDEPEND) $(CFLAGS)  $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) $(LOCAL_INCLUDES) *.cc > make.depend

include make.depend				
