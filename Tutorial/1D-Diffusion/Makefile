# This is the Loci base directory.
LOCI_BASE ?= /opt/loci
# List the name of your compiled program here
TARGET = heat

include $(LOCI_BASE)/Loci.conf
#COPT=-g
# List your object files here
OBJS  =  main.o heat.o
OBJS2 =  main.o heat2.o

JUNK = *~  core make.depend make.depend.bak

default:
	$(MAKE) $(TARGET)

$(TARGET): $(OBJS) FRC
	$(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

heat2: $(OBJS2) FRC
	$(LD) -o heat2 $(OBJS2) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

FRC : 

clean:
	rm -fr $(OBJS) $(OBJS2) $(TARGET) heat2 $(JUNK) 

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)
DEPEND_FILES=$(subst .o,.d,$(OBJS))

distclean:
	rm -fr $(OBJS) $(OBJS2) $(TARGET) heat2 $(JUNK) $(LOCI_LPP_FILES) $(DEPEND_FILES) debug



#include automatically generated dependencies
include $(DEPEND_FILES)
