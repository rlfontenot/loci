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

DIR_NAME	=	FVMMod
MODULE_NAME	=	fvm
LIB_NAME	=	$(MODULE_NAME)_m.so

################################################################################
# No changes required below this line
################################################################################
CXX_FILES		=	$(wildcard *.cpp)
LOCI_FILES		=	$(wildcard *.loci)
LOCI_LPP_FILES	=	$(LOCI_FILES:.loci=.cc)
OBJS			=	$(subst .loci,.o,$(LOCI_FILES)) \
					$(subst .cpp,.o,$(CXX_FILES))
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))
JUNK			=	*~ core ti_files ii_files rii_files
DEPS			=	$(LOCI_SRC)/bin/lpp \
					$(LOCI_SRC)/lib/libLoci.so
INCLUDES		=	-I$(LOCI_SRC)/include


default:	$(LIB_NAME)
all:		$(LIB_NAME)
spotless:	distclean clean_install

$(LOCI_LPP_FILES):	$(DEPS)

$(LOCI_SRC)/bin/lpp:
	$(MAKE) -C $(LOCI_SRC)/lpp
$(LOCI_SRC)/lib/libLoci.so:
	$(MAKE) -C $(LOCI_SRC)/System

$(LIB_NAME): $(OBJS)
	@echo -e "\e[$(PRT) $@ \e[0m"
	$(SHARED_LD) $(SHARED_LD_FLAGS) $@ $(LDFLAGS) $(LDLIBS) $(OBJS)
	cp $@ $(LOCI_SRC)/lib/$@

install: $(LOCI_BASE)/include/$(DIR_NAME)/. $(LOCI_BASE)/lib/. $(LIB_NAME)
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	cp $(LOCI_SRC)/lib/$(LIB_NAME) $(LOCI_BASE)/lib/$(LIB_NAME)
	cp -u $(LOCI_SRC)/include/$(DIR_NAME)/*.h  $(LOCI_BASE)/include/$(DIR_NAME)/
	cp -u $(LOCI_SRC)/include/$(DIR_NAME)/*.lh $(LOCI_BASE)/include/$(DIR_NAME)/
endif

clean_install:
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME)
	rm -rf $(LOCI_BASE)/include/$(DIR_NAME)
endif

clean:
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"
	rm -fr $(LIB_NAME) $(OBJS) $(LOCI_SRC)/lib/$(LIB_NAME) $(JUNK)

distclean: clean
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"
	rm -fr  $(LOCI_LPP_FILES) $(DEPEND_FILES)

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif
