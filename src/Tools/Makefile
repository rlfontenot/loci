#
# Tools Makefile
# 

TARGET = libTools.a

LOCAL_INCLUDES = -I$(LOCI_BASE)/include

OBJS = fpe.o ftrn_reader.o identity.o eventNotify.o parse.o options_list.o \
       expr.o intervalSet.o digraph.o


.SUFFIXES: .o .c .cc .f

.c.o:    ;  $(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.c
.cc.o:   ;  $(CP) $(CFLAGS) $(TEMPLATE_LIB_OPT) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.cc
.f.o:    ;  $(F77) $(FFLAGS) -c $*.f

default:
	cd .. ; make -f Makefile Tools

compile: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS) 
	$(RANLIB) $(TARGET)



JUNK =  *~ core ti_files KCC_files so_locations

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)


make.depend:
	touch make.depend 

depend: make.depend
	makedepend -fmake.depend $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -I$(CC_INCLUDE_PATH) *.cc

include make.depend				
