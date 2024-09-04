################################################################################
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
################################################################################
LOCI_BASE?=$(shell pwd)
LOCI_SRC?=$(shell pwd)

include $(LOCI_SRC)/Loci.conf
include $(LOCI_SRC)/targets/external.mk

MODELS_BIN	=	2dgv          \
				ccm2vog       \
				ccm2vogp      \
				cfd++2vog     \
				cgns2ensight  \
				cgns2surf     \
				cgns2ugrid    \
				cgns2vog      \
				cobalt2vog    \
				extract       \
				extract_chem  \
				extruder      \
				fluent2vog    \
				lpp           \
				make_periodic \
				marker        \
				msh2vog       \
				plot3d2vog    \
				refine        \
				refmesh       \
				ugrid2cgns    \
				ugrid2vog     \
				ugridmirror   \
				vog2surf      \
				vogmerge      \
				vogcut        \
				vogcheck

MODELS_LIB	=	fvmadapt_m.so                 \
				fvmadaptFAD_m.so              \
				fvm_m.so                      \
				fvmFAD_m.so                   \
				fvmoverset_m.so               \
				fvmoversetFAD_m.so            \
				libadf.a                      \
				libadf.$(LIB_SUFFIX)          \
				libfvmadaptfunc.$(LIB_SUFFIX) \
				libLoci.$(LIB_SUFFIX)         \
				libsprng.$(LIB_SUFFIX)        \
				libTools.a                    \
				libTools.$(LIB_SUFFIX)

default:
	-$(MAKE) -C 2dgv         $@
	$(MAKE)  -C Tools        $@
	$(MAKE)  -C System       $@
	$(MAKE)  -C lpp          $@
	$(MAKE)  -C FVMAdapt     $@
	$(MAKE)  -C FVMMod       $@
	$(MAKE)  -C FVMOverset   $@
	$(MAKE)  -C FVMtools     $@
	$(MAKE)  -C sprng        $@
	$(MAKE)  -C Tutorial     $@

copy_config:
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	mkdir -p $(LOCI_BASE)
	$(COPY) *.conf $(LOCI_BASE)/
	$(SED_I) "s,LOCI_SRC,LOCI_BASE,g" $(LOCI_BASE)/Loci.conf
	$(SED_I) "s,LOCI_SRC,LOCI_BASE,g" $(LOCI_BASE)/sys.conf
endif

install: copy_config
	$(MAKE)  $@_ParMETIS
	$(MAKE)  $@_METIS
	$(MAKE)  $@_GKlib
	-$(MAKE) -C 2dgv         $@
	$(MAKE)  -C FVMAdapt     $@
	$(MAKE)  -C FVMMod       $@
	$(MAKE)  -C FVMOverset   $@
	$(MAKE)  -C FVMtools     $@
	$(MAKE)  -C lpp          $@
	$(MAKE)  -C sprng        $@
	$(MAKE)  -C System       $@
	$(MAKE)  -C Tools        $@
	$(MAKE)  -C Tutorial     $@

clean:
	$(MAKE) -C 2dgv         $@
	$(MAKE) -C FVMAdapt     $@
	$(MAKE) -C FVMMod       $@
	$(MAKE) -C FVMOverset   $@
	$(MAKE) -C FVMtools     $@
	$(MAKE) -C lpp          $@
	$(MAKE) $@_ParMETIS
	$(MAKE) $@_METIS
	$(MAKE) $@_GKlib
	$(MAKE) -C sprng        $@
	$(MAKE) -C System       $@
	$(MAKE) -C Tools        $@
	$(MAKE) -C Tutorial     $@

distclean:
	$(MAKE) -C 2dgv         $@
	$(MAKE) -C FVMAdapt     $@
	$(MAKE) -C FVMMod       $@
	$(MAKE) -C FVMOverset   $@
	$(MAKE) -C FVMtools     $@
	$(MAKE) -C lpp          $@
	$(MAKE) -C sprng        $@
	$(MAKE) -C System       $@
	$(MAKE) -C Tools        $@
	$(MAKE) -C Tutorial     $@

clean_install:
	$(MAKE) -C 2dgv         $@
	$(MAKE) -C FVMAdapt     $@
	$(MAKE) -C FVMMod       $@
	$(MAKE) -C FVMOverset   $@
	$(MAKE) -C FVMtools     $@
	$(MAKE) -C lpp          $@
	$(MAKE) -C sprng        $@
	$(MAKE) -C System       $@
	$(MAKE) -C Tools        $@
	$(MAKE) -C Tutorial     $@
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/Loci.conf
	rm -rf $(LOCI_BASE)/sys.conf
	rm -rf $(LOCI_BASE)/comp.conf
endif

spotless:
	$(MAKE) -C 2dgv         $@
	$(MAKE) -C FVMAdapt     $@
	$(MAKE) -C FVMMod       $@
	$(MAKE) -C FVMOverset   $@
	$(MAKE) -C FVMtools     $@
	$(MAKE) -C lpp          $@
	$(MAKE) -C sprng        $@
	$(MAKE) -C System       $@
	$(MAKE) -C Tools        $@
	$(MAKE) -C Tutorial     $@
	$(MAKE) $@_ParMETIS
	$(MAKE) $@_METIS
	$(MAKE) $@_GKlib
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/Loci.conf
	rm -rf $(LOCI_BASE)/sys.conf
	rm -rf $(LOCI_BASE)/comp.conf
endif


RED		=	$(ESC)1;31;1m
GREEN	=	$(ESC)1;32;1m
DIREC	=	$(LOCI_SRC)
list:
	@$(ECHO) "$(PRT) :: $(DIREC) :: $(NC)"
	@$(ECHO) "$(PRT) :: bin :: $(NC)"
	@for d in $(MODELS_BIN); do \
		if [ -f ${DIREC}/bin/$$d ]; then \
			$(ECHO) "$(GREEN) $$d $(NC)"; \
		else \
			$(ECHO) "$(RED) $$d $(NC)"; \
		fi \
	done
	@$(ECHO) "$(PRT) :: lib :: $(NC)"
	@for d in $(MODELS_LIB); do \
		if [ -f ${DIREC}/lib/$$d ]; then \
			$(ECHO) "$(GREEN) $$d $(NC)"; \
		else \
			$(ECHO) "$(RED) $$d $(NC)"; \
		fi \
	done

list_installed: DIREC=$(LOCI_BASE)
list_installed: list


PRTP	=	$(ESC)1;35;1m
devhelp:
	@$(ECHO) "$(PRTP)LOCI_SRC           $(NC) = $(LOCI_SRC)"
	@$(ECHO) "$(PRTP)LOCI_BASE          $(NC) = $(LOCI_BASE)"
	@$(ECHO) "$(PRTP)CC                 $(NC) = $(CC)"
	@$(ECHO) "$(PRTP)CXX                $(NC) = $(CXX)"
	@$(ECHO) "$(PRTP)FC                 $(NC) = $(FC)"
	@$(ECHO) "$(PRTP)MAKEDEPEND         $(NC) = $(MAKEDEPEND)"
	@$(ECHO) "$(PRTP)LD                 $(NC) = $(LD)"
	@$(ECHO) "$(PRTP)SHARED_LD          $(NC) = $(SHARED_LD)"
	@$(ECHO) "$(PRTP)DYNAMIC_LD         $(NC) = $(DYNAMIC_LD)"
	@$(ECHO) "$(PRTP)FFLAGS             $(NC) = $(FFLAGS)"
	@$(ECHO) "$(PRTP)CFLAGS             $(NC) = $(CFLAGS)"
	@$(ECHO) "$(PRTP)CXXFLAGS           $(NC) = $(CXXFLAGS)"
	@$(ECHO) "$(PRTP)LDFLAGS            $(NC) = $(LDFLAGS)"
	@$(ECHO) "$(PRTP)SHARED_LD_FLAGS    $(NC) = $(SHARED_LD_FLAGS)"
	@$(ECHO) "$(PRTP)DYNAMIC_LD_FLAGS   $(NC) = $(DYNAMIC_LD_FLAGS)"
	@$(ECHO) "$(PRTP)LDPATH             $(NC) = $(LDPATH)"
	@$(ECHO) "$(PRTP)RPATH              $(NC) = $(RPATH)"
	@$(ECHO) "$(PRTP)DEFINES            $(NC) = $(DEFINES)"
	@$(ECHO) "$(PRTP)FAD_DEFINES        $(NC) = $(FAD_DEFINES)"
	@$(ECHO) "$(PRTP)LOCI_INCLUDES      $(NC) = $(LOCI_INCLUDES)"
	@$(ECHO) "$(PRTP)INCLUDES           $(NC) = $(INCLUDES)"
	@$(ECHO) "$(PRTP)LDLIBS             $(NC) = $(LDLIBS)"
	@$(ECHO) "$(PRTP)LOCAL_LIBS         $(NC) = $(LOCAL_LIBS)"
	@$(ECHO) "$(PRTP)LOCI_DEBUG         $(NC) = $(LOCI_DEBUG)"


