# This is the Loci base directory.  The one listed here defines the
# debugging version for the sgi platform

LOCI_BASE = /ccs/mssl/general/loci/opt/Loci-sgi-rel-0-3-4f

# List the name of your compiled program here
TARGET = fire

# List your object files here
OBJS =  main.o advance_time.o boundary_conditions.o cell_center.o \
        collapse_time.o gradient.o initial_conditions.o stable_timestep.o


include $(LOCI_BASE)/sys.conf
include $(LOCI_BASE)/comp.conf

DEFINES = $(SYSTEM_DEFINES) $(MACHINE_SPECIFIC) $(DEBUG) -DBOUNDS_CHECK

INCLUDES  = -I$(LOCI_BASE)/include 
LIBS      = -L$(LOCI_BASE)/lib -lLoci -lTools $(BASE_LIBS) 

CP  = $(CPP) $(CPPOPT) -w 
CC  = $(CPP) $(CPPOPT) -w 

COPT     = $(CC_DEBUG) $(CC_OPT)
COPTLESS = $(CC_DEBUG) $(CC_OPT1)

LOPT     = $(CC_DEBUG) $(COPT)

JUNK = *~  core ti_files ii_files make.depend make.depend.bak

.SUFFIXES: .o .c .cc

.c.o:    ;  $(CC) $(COPT) $(DEFINES) $(INCLUDES) -c $*.c
.cc.o:   ;  $(CP) $(COPT) $(DEFINES) $(INCLUDES) -c $*.cc

default:
	$(MAKE) $(TARGET)

$(TARGET): $(OBJS) FRC
	$(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

FRC : 

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)  rii_files ii_files

