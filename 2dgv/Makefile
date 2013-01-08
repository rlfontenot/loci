# Put the object filenames here 
# (replace AUTOMATIC_OBJS with list of .o files if you don't want to compile
# all files in this directory into a single executable)
OBJS = $(AUTOMATIC_OBJS) 

#Base directory where X11 include files and libraries are located
X11BASE=/usr/X11R6

# Put the executable name here
TARGET = 2dgv
# Put linker flags here (such as any libraries to link)
LIBRARIES = -lm -L$(X11BASE)/lib -lXt -lXaw -lX11 

# Include search directories
INCLUDES = -I. -I$(X11BASE)/include

# Compiler flags
# C Compiler
#CC = gcc
# Put C compiler flags here (default debugging options, basic optimization)
CFLAGS= -O1 $(INCLUDES)

# C++ Compiler
#CXX = g++
# Put C++ Compiler Flags here (default debugging options, basic optimization)
CXXFLAGS= -O1 $(INCLUDES)



#############################################################################
# No need to change rules below this line
#############################################################################

# Find program files in this directory
AUTOMATIC_FILES = $(wildcard *.c *.cc *.C)
AUTOMATIC_OBJS = $(subst .c,.o,$(subst .cc,.o,$(subst .C,.o,$(AUTOMATIC_FILES))))

# Compile target program
$(TARGET): $(OBJS) 
	$(CXX) -o $(TARGET) $(OBJS) $(LIBRARIES)


# rule for generating dependencies from source files
%.d: %.c
	set -e; $(CC) -M $(CFLAGS) $< \
	| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
		[ -s $@ ] || rm -f $@
%.d: %.C
	set -e; $(CXX) -M $(CXXFLAGS) $< \
	| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
		[ -s $@ ] || rm -f $@
%.d: %.cc
	set -e; $(CXX) -M $(CXXFLAGS) $< \
	| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
		[ -s $@ ] || rm -f $@

DEPEND_FILES=$(subst .o,.d,$(OBJS))


clean:
	rm -f $(OBJS) $(TARGET) 

distclean:
	rm -f $(OBJS) $(TARGET) $(DEPEND_FILES)

#include automatically generated dependencies
include $(DEPEND_FILES)


