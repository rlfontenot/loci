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

DIR_NAME	=	FVMOverset
MODULE_NAME	=	fvmoverset
LIB_NAME	=	$(MODULE_NAME)_m.so

################################################################################
# No changes required below this line
################################################################################
LOCI_FILES		=	$(wildcard *.loci)
LOCI_LPP_FILES	=	$(LOCI_FILES:.loci=.cc)
OBJS			=	$(subst .loci,.o,$(LOCI_FILES))
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))
JUNK			=	*~  core ti_files ii_files rii_files
INCLUDES		=	-I$(LOCI_SRC)/include/FVMOverset \
					-I$(LOCI_SRC)/include
DEPS			=	$(LOCI_SRC)/bin/lpp


default:	$(LIB_NAME)
all:		$(LIB_NAME)
spotless:	distclean clean_install

$(LOCI_LPP_FILES): $(DEPS)

$(LOCI_SRC)/bin/lpp:
	$(MAKE) -C $(LOCI_SRC)/lpp


$(LIB_NAME): $(OBJS)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(SHARED_LD) $(SHARED_LD_FLAGS) $@ $(LDFLAGS) $(OBJS)
	$(COPY) $@ $(LOCI_SRC)/lib/$@

install: $(LOCI_BASE)/lib/. $(LOCI_BASE)/include/$(DIR_NAME)/. $(LIB_NAME)
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	$(COPY) $(LOCI_SRC)/lib/$(LIB_NAME) $(LOCI_BASE)/lib/$(LIB_NAME)
	$(COPY) $(LOCI_SRC)/include/$(DIR_NAME)/*.h $(LOCI_BASE)/include/$(DIR_NAME)/
	$(COPY) $(LOCI_SRC)/include/$(DIR_NAME)/*.lh $(LOCI_BASE)/include/$(DIR_NAME)/
endif

clean_install:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME)
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/include/$(DIR_NAME)
endif

clean:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(OBJS) $(LIB_NAME) $(LOCI_SRC)/lib/$(LIB_NAME) $(JUNK)

distclean: clean
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr  $(LOCI_LPP_FILES) $(DEPEND_FILES)


#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif
