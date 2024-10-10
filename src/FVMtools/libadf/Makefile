
include $(LOCI_BASE)/Loci.conf

TARGET = libadf.a

INCLUDES=-I.

COPT = $(CC_LIB_FLAGS) $(COPTLESS)
C_OPT = $(CC_LIB_FLAGS) $(COPTLESS)

OBJS = ADF_interface.o ADF_internals.o

all: $(TARGET)

libadf.a: $(OBJS) 
	$(AR) $(ARFLAGS) libadf.a $(OBJS) 
	$(RANLIB) libadf.a

LIB_OBJS=$(OBJS:.o=_lo.o)
libadf.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) libadf.so $(LIB_FLAGS) $(LIB_OBJS) 

JUNK =  *~ core ii_files ti_files so_locations rii_files

clean:
	rm -fr $(LIB_OBJS) $(OBJS) $(TARGET) $(JUNK)


distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) #$(DEPEND_FILES)
	

#DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
#-include $(DEPEND_FILES)
# DO NOT DELETE
