# LOCI_BASE should be set before the Makefile

# List your object files here
OBJS  =  bcFixedMass.o bcImpermeable.o bcOutflow.o gridReader.o gridSetup.o gridTools.o fluxInviscid.o main.o IOgrid.o IOplot.o solverPeriodic.o solverRungeKutta.o solverSetup.o solverTimestep.o solverIC.o solverEoS.o
INCLUDES = -I./include

# List the name of your compiled program here
TARGET = fluidTest


#############################################################################
# No changes needed below this line

include $(LOCI_BASE)/Loci.conf
#COPT= -g -O0

default: fluidTest ndiff
	cd tests; $(MAKE) TestResults; cd ..

fluidTest: $(OBJS)
	$(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

ndiff: ndiff.o
	$(LD) -o  ndiff ndiff.o $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

clean:
	rm -fr $(OBJS) $(TARGET)
	cd tests; $(MAKE) clean; cd ..

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~) $(wildcard crash_dump.*)  core debug output $(wildcard *.o) ndiff
# ".cc" files created from .loci files
LOCI_INTERMEDIATE_FILES = $(subst .loci,.cc, $(wildcard *.loci) )

distclean:
	rm -fr $(OBJS) $(TARGET) $(JUNK) $(LOCI_INTERMEDIATE_FILES) $(DEPEND_FILES)
	cd tests; $(MAKE) distclean; cd ..

DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
-include $(DEPEND_FILES)
