#
# System Makefile
# 

TARGET = libLoci.a #libLoci.so

LOCAL_INCLUDES = -I$(LOCI_BASE)/include


OBJS = variable.o store_rep.o Map.o constraint.o rule.o \
       fact_db.o execute.o depend_graph.o scheduler.o


JUNK = *~ *.ii core ti_files ii_files 

.SUFFIXES: .o .c .cc .f 

.c.o:    ;  $(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.c
.cc.o:   ;  $(CP) $(CFLAGS) $(TEMPLATE_LIB_OPT) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c $*.cc
.f.o:    ;  $(F77) $(FFLAGS) -c $*.f


default:
	cd .. ; make -f Makefile System


compile: $(TARGET)

libLoci.a: $(OBJS)
	$(AR) $(ARFLAGS) libLoci.a $(OBJS) 
	$(RANLIB) $(TARGET)

libLoci.so: $(OBJS)
	$(LD) $(LDFLAGS) -shared -o libLoci.so $(OBJS)

scheduler.o: scheduler.cc
	$(CP) $(CFLAGS) $(TEMPLATE_LIB_OPT) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES) -c scheduler.cc

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)


make.depend:
	touch make.depend 

depend: make.depend
	$(MAKEDEPEND) $(DEFINES) $(INCLUDES) $(LOCAL_INCLUDES)  *.cc > make.depend

include make.depend				
