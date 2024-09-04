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

DIR_NAME	=	System
MODULE_NAME	=	Loci
LIB_NAME1	=	lib$(MODULE_NAME).a
LIB_NAME2	=	lib$(MODULE_NAME).$(LIB_SUFFIX)

################################################################################
# No changes required below this line
################################################################################
# Remove the filter when gpu stuff is fully developed
FILES			=	$(filter-out \
					gconstraint.cpp \
					gkey_manager.cpp\
					gkeyspace.cpp \
					gmap.cpp \
					gmultimap.cpp \
					,$(wildcard *.cpp))
OBJS			=	$(subst .cpp,.o,$(FILES))
LS_INC_H		=	$(shell cd $(LOCI_SRC)/include/; ls *.h) version_loci.h
LS_INC_LH		=	$(shell cd $(LOCI_SRC)/include/; ls *.lh)
GIT_INFO		=	"$(shell git describe --abbrev=6 --always --dirty --tags)"
GIT_BRANCH		=	"$(shell git rev-parse --abbrev-ref HEAD)"
VERSION_FILE	=	${LOCI_SRC}/include/version_loci.h
DATE			=	$(shell date +"%A %b. %d %Y, %r %Z")
OLD_GIT_INFO	=	"$(shell cat $(VERSION_FILE) 2>/dev/null | \
                       grep LOCI_VERSION                     | \
                       cut -d'"' -f2)"
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))
JUNK			=	*~ core ii_files ti_files KCC_files so_locations rii_files
DEFINES			+=	-DLOCI_RPATH=\"$(LOCI_RPATH)\"
INCLUDES		=	-I$(LOCI_SRC)/include/sprng \
					-I$(LOCI_SRC)/include       \
					$(PARMETIS_INCLUDE)
LOCAL_LIBS		=	-lparmetis -lsprng
DEPS			=	${PARMETIS_BASE}/lib/libparmetis.$(LIB_SUFFIX) \
					${LOCI_SRC}/lib/libsprng.$(LIB_SUFFIX)


default:	$(LIB_NAME2)
all:		$(LIB_NAME1) $(LIB_NAME2)
spotless:	distclean clean_install

${PARMETIS_BASE}/lib/libparmetis.$(LIB_SUFFIX):
	@$(MAKE) -C ${LOCI_SRC} ParMETIS
${LOCI_SRC}/lib/libsprng.$(LIB_SUFFIX):
	@$(MAKE) -C ${LOCI_SRC}/sprng all

FVMGridReader.o: $(DEPS)

$(LIB_NAME1): $(DEPS) $(VERSION_FILE) $(OBJS)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(RANLIB) $@
	$(COPY) -u $@ $(LOCI_SRC)/lib/$@

$(LIB_NAME2): $(DEPS) $(VERSION_FILE) $(OBJS)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(DYNAMIC_LD) $(DYNAMIC_LD_FLAGS) $@ $(LDFLAGS) $(LOCAL_LIBS) $(OBJS)
	$(COPY) -u $@ $(LOCI_SRC)/lib/$@

install: $(LIB_NAME2) install_include
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	$(COPY) -u $(LOCI_SRC)/lib/$(LIB_NAME2) $(LOCI_BASE)/lib/$(LIB_NAME2)
endif

install_include: $(LOCI_BASE)/include/Config/.
install_include: $(LOCI_BASE)/include/MPI_stubb/.
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	$(COPY) -u $(LOCI_SRC)/include/Config/*.h    $(LOCI_BASE)/include/Config/
	$(COPY) -u $(LOCI_SRC)/include/MPI_stubb/*.h $(LOCI_BASE)/include/MPI_stubb/
	$(COPY) -u $(LOCI_SRC)/include/Loci $(LOCI_BASE)/include/
	for d in $(LS_INC_H);  do $(COPY) -p  $(LOCI_SRC)/include/$$d $(LOCI_BASE)/include/; done
	for d in $(LS_INC_LH); do $(COPY) -p  $(LOCI_SRC)/include/$$d $(LOCI_BASE)/include/; done
endif

clean_install: clean_install_include
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME1)
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME2)

clean_install_include:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	for d in $(LS_INC_H);  do rm -f $(LOCI_BASE)/include/$$d; done
	for d in $(LS_INC_LH); do rm -f $(LOCI_BASE)/include/$$d; done
	rm -rf $(LOCI_BASE)/include/Config
	rm -rf $(LOCI_BASE)/include/MPI_stubb
	rm -rf $(LOCI_BASE)/include/Loci
endif

clean:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(OBJS) $(LIB_NAME1) $(LIB_NAME2) $(JUNK)
	rm -fr $(LOCI_SRC)/lib/$(LIB_NAME1)
	rm -fr $(LOCI_SRC)/lib/$(LIB_NAME2)

distclean: clean
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(DEPEND_FILES)


$(VERSION_FILE):
ifneq ($(shell date '+%F %H:%M' -r $@ 2>/dev/null),$(shell date \
                                                                  '+%F %H:%M'))
ifneq ($(GIT_INFO),$(OLD_GIT_INFO))
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	@rm -rf $@
	@echo "/** ****************************************************************\
	************"                                    >> $@
	@echo " * @file    version_loci.h"               >> $@
	@echo " * @author  Mark A. Hunt (CFDRC)"         >> $@
	@echo " * @brief   Auto-generated version file." >> $@
	@echo " * @version 0.1"                          >> $@
	@echo " * @date    2024-06-18"                   >> $@
	@echo " * @copyright Copyright (c) 2024"         >> $@
	@echo " *******************************************************************\
	***********/"                                    >> $@
	@echo ""                                         >> $@
	@echo "#ifndef VERSION_LOCI_H"                   >> $@
	@echo "#define VERSION_LOCI_H"                   >> $@
	@echo ""                                         >> $@
	@echo "#define LOCI_BRANCH   \"$(GIT_BRANCH)\""  >> $@
	@echo "#define LOCI_VERSION  \"$(GIT_INFO)\""    >> $@
	@echo "#define LOCI_DATETIME \"$(DATE)\""        >> $@
	@echo ""                                         >> $@
	@echo "#endif"                                   >> $@
endif
endif

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif
