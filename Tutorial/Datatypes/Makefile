#LOCI_BASE = /ccs/mssl/general/loci/debug/sparc/
LOCI_BASE = /usr/tmp/csv/newsrc/Loci
TARGET = fire

# List your object files here

OBJS = testDMultiMap.o

include $(LOCI_BASE)/sys.conf
include $(LOCI_BASE)/comp.conf

DEFINES = $(SYSTEM_DEFINES) $(MACHINE_SPECIFIC) $(DEBUG) -DBOUNDS_CHECK

INCLUDES  = -I$(LOCI_BASE)/include 
LIBS      = -L/usr/local/lib -lz -L$(LOCI_BASE)/lib -lLoci -lTools -lhdf5CC $(BASE_LIBS) \
            -L/ccs/mssl/general/loci/misc/sparc/lib -lmetis -lm 

CP  = $(CPP) $(CPPOPT) 
CC  = $(CPP) $(CPPOPT)

COPT     = $(CC_DEBUG) $(CC_OPT)
COPTLESS = $(CC_DEBUG) $(CC_OPT1)

JUNK = *~ *.ii core q.file_1.* ti_files 

.SUFFIXES: .o .c .cc

.c.o:    ;  $(CC) -w $(CC_OPT) $(DEFINES) $(INCLUDES) -c $*.c
.cc.o:   ;  $(CP) -w $(CC_OPT) $(DEFINES) $(INCLUDES) -c $*.cc

default:
	$(MAKE) $(TARGET)

$(TARGET): $(OBJS) FRC
	$(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

newgrid: newgrid.o
	$(LD) -o newgrid newgrid.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

FRC : 

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK) make.depend make.depend.bak

