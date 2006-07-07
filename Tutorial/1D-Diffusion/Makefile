# This is the Loci base directory.
LOCI_BASE ?= /usr/tmp/lush/Loci
# List the name of your compiled program here
TARGET = onedheat

include $(LOCI_BASE)/Loci.conf
#COPT=-g
# List your object files here
OBJS =  main.o heat2.o


JUNK = *~  core ti_files ii_files make.depend make.depend.bak

default:
	$(MAKE) $(TARGET)

$(TARGET): $(OBJS) FRC
	$(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

FRC : 

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)  rii_files ii_files

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
include $(DEPEND_FILES)
