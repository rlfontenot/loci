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
LOCI_BASE?=..
LOCI_SRC?=..

include $(LOCI_SRC)/Loci.conf
TARGET	=	lpp

################################################################################
# No changes required below this line
################################################################################
FILES			=	$(wildcard *.cpp)
OBJS			=	$(subst .cpp,.o,$(FILES)) \
					../System/variable.o
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))
JUNK			=	*~ core ii_files ti_files KCC_files so_locations rii_files
INCLUDES		=	-I$(LOCI_SRC)/include
LOCAL_LIBS		=	$(RPATH)$(LOCI_RPATH) -L$(LOCI_SRC)/lib -lTools -lm
DEPS			=	$(LOCI_SRC)/lib/libTools.so

default:	$(TARGET)
all:		$(TARGET)
spotless:	distclean clean_install


$(LOCI_SRC)/lib/libTools.so:
	$(MAKE) -C $(LOCI_SRC)/Tools

$(TARGET): $(DEPS) $(OBJS)
	@echo -e "\e[$(PRT) $@ \e[0m"
	$(LDBSL) -o $@ $(OBJS) $(LOCAL_LIBS) $(LDFLAGS)
	cp $@ $(LOCI_SRC)/bin/$@

install: $(LOCI_BASE)/bin/. $(TARGET)
	@echo -e "\e[$(PRT) $@_$(TARGET) \e[0m"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	cp -u $(LOCI_SRC)/bin/$(TARGET) $(LOCI_BASE)/bin/$(TARGET)
endif

clean_install:
	@echo -e "\e[$(PRT) $@_$(TARGET) \e[0m"
	rm -rf $(LOCI_BASE)/bin/$(TARGET)

clean:
	@echo -e "\e[$(PRT) $@_$(TARGET) \e[0m"
	rm -fr $(OBJS) $(TARGET) $(LOCI_SRC)/bin/$(TARGET) $(JUNK)

distclean: clean
	@echo -e "\e[$(PRT) $@_$(TARGET) \e[0m"
	rm -fr $(DEPEND_FILES)


ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif
