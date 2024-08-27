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

MODULE_NAME	=	Tools
LIB_NAME1	=	lib$(MODULE_NAME).a
LIB_NAME2	=	lib$(MODULE_NAME).$(LIB_SUFFIX)

################################################################################
# No changes required below this line
################################################################################
CFILES			=	$(wildcard *.c)
CXXFILES		=	$(wildcard *.cpp)
OBJS			=	$(subst .c,.o,$(CFILES)) \
					$(subst .cpp,.o,$(CXXFILES))
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))
JUNK			=	*~ core ii_files ti_files KCC_files so_locations rii_files
INCLUDES		=	-I$(LOCI_SRC)/include/Tools \
					-I$(LOCI_SRC)/include


default:	$(LIB_NAME1) $(LIB_NAME2)
all:		$(LIB_NAME1) $(LIB_NAME2)
spotless:	distclean clean_install

$(LIB_NAME1): $(OBJS)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(RANLIB) $@
	$(COPY) -u $@ $(LOCI_SRC)/lib/$@

$(LIB_NAME2): $(OBJS)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(DYNAMIC_LD) $(DYNAMIC_LD_FLAGS) $@ $(LDFLAGS) $(OBJS)
	$(COPY) -u $@ $(LOCI_SRC)/lib/$@


install: $(LOCI_BASE)/lib/. $(LOCI_BASE)/include/$(MODULE_NAME)/. $(LIB_NAME1) $(LIB_NAME2)
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	$(COPY) -u $(LOCI_SRC)/lib/$(LIB_NAME1)        $(LOCI_BASE)/lib/$(LIB_NAME1)
	$(COPY) -u $(LOCI_SRC)/lib/$(LIB_NAME2)        $(LOCI_BASE)/lib/$(LIB_NAME2)
	$(COPY) -u $(LOCI_SRC)/include/$(MODULE_NAME)/*.h $(LOCI_BASE)/include/$(MODULE_NAME)/
endif

clean_install:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME1)
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME2)
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/include/$(MODULE_NAME)
endif

clean:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(OBJS) $(LIB_NAME1) $(LIB_NAME2)
	rm -fr $(LOCI_SRC)/lib/$(LIB_NAME1) $(LOCI_SRC)/lib/$(LIB_NAME2) $(JUNK)

distclean: clean
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(DEPEND_FILES)

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif

