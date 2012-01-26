###############################################################################
#
# Copyright 2008, Mississippi State University
#
# This file is part of the Loci Framework.
#
# The Loci Framework is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The Loci Framework is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
#
###############################################################################
include $(LOCI_BASE)/Loci.conf
#COPT = -g

# List the name of your compiled program here
TARGET = lpp


# List your object files here
OBJS =  lpp.o main.o variable.o


LOCAL_LIBS =

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~)  core 

default:
	$(MAKE) $(TARGET) 

all: $(TARGET) 


TOOLS_LIBS = $(RPATH)$(LOCI_RPATH) -L$(LOCI_BASE)/lib -lTools -lm

$(TARGET): $(OBJS) 
	$(LDBSL) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(TOOLS_LIBS) $(LDFLAGS)
FRC : 

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK)

distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) $(DEPEND_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
