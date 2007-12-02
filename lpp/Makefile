#LOCI_BASE=/usr/local/Loci-Linux-x86_64-mpic++-rel-3-0-p8/
include $(LOCI_BASE)/Loci.conf
#COPT = -g

# List the name of your compiled program here
TARGET = lpp


# List your object files here
OBJS =  lpp.o main.o


LOCAL_LIBS =

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~)  core 

default:
	$(MAKE) $(TARGET) 

all: $(TARGET) 


$(TARGET): $(OBJS) 
	$(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)
FRC : 

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK)

distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
