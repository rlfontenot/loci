#
# System Makefile
# 

TARGET = libLoci.a

LOCAL_INCLUDES = -I$(LOCI_BASE)/include -I$(LOCI_BASE)/include/Tools


OBJS = variable.o store_rep.o Map.o constraint.o path.o rule.o \
       fact_db.o execute.o depend_graph.o scheduler.o


JUNK = *~ *.ii core ti_files ii_files Templates.DB KCC_files so_locations

.SUFFIXES: .o .c .cc .f 

.c.o:    ;  $(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.c
.cc.o:   ;  $(CP) $(CFLAGS) $(TEMPLATE_LIB_OPT) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.cc
.f.o:    ;  $(F77) $(FFLAGS) -c $*.f


default:
	cd .. ; make -f Makefile System


compile: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS) 
	$(RANLIB) $(TARGET)

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)


make.depend:
	touch make.depend 

depend: make.depend
	makedepend -fmake.depend $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -I$(CC_INCLUDE_PATH) *.cc

include make.depend				
