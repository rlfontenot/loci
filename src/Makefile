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

MODELS_LIB	=	fvmadapt_m.so      \
				fvm_m.so           \
				fvmoverset_m.so    \
				libadf.a           \
				libadf.so          \
				libfvmadaptfunc.so \
				libLoci.so         \
				libsprng.so        \
				libTools.a         \
				libTools.so

default:
	$(MAKE) -C 2dgv         $@
	$(MAKE) -C Tools        $@
	$(MAKE) -C System       $@
	$(MAKE) -C lpp          $@
	$(MAKE) -C FVMAdapt     $@
	$(MAKE) -C FVMMod       $@
	$(MAKE) -C FVMOverset   $@
	$(MAKE) -C FVMtools     $@
	$(MAKE) -C sprng        $@
	$(MAKE) -C Tutorial     $@

copy_config:
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	mkdir -p $(LOCI_BASE)
	cp *.conf $(LOCI_BASE)/
	sed -i "s,LOCI_SRC,LOCI_BASE,g" $(LOCI_BASE)/Loci.conf
	sed -i "s,LOCI_SRC,LOCI_BASE,g" $(LOCI_BASE)/sys.conf
endif

install: copy_config
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
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/Loci.conf
	rm -rf $(LOCI_BASE)/sys.conf
	rm -rf $(LOCI_BASE)/comp.conf
endif


RED		=	\e[1;31;1m
GREEN	=	\e[1;32;1m
DIREC	=	$(LOCI_SRC)
list:
	@echo -e "\e[$(PRT) :: $(DIREC) :: \e[0m"
	@echo -e "\e[$(PRT) :: bin :: \e[0m"
	@for d in $(MODELS_BIN); do \
		if [ -f ${DIREC}/bin/$$d ]; then \
			echo -e "$(GREEN) $$d \e[0m"; \
		else \
			echo -e "$(RED) $$d \e[0m"; \
		fi \
	done
	@echo -e "\e[$(PRT) :: lib :: \e[0m"
	@for d in $(MODELS_LIB); do \
		if [ -f ${DIREC}/lib/$$d ]; then \
			echo -e "$(GREEN) $$d \e[0m"; \
		else \
			echo -e "$(RED) $$d \e[0m"; \
		fi \
	done

list_installed: DIREC=$(LOCI_BASE)
list_installed: list


PRTP	=	\e[1;35;1m
devhelp:
	@echo -e "$(PRTP)LOCI_SRC           \e[0m = $(LOCI_SRC)"
	@echo -e "$(PRTP)LOCI_BASE          \e[0m = $(LOCI_BASE)"
	@echo -e "$(PRTP)CC                 \e[0m = $(CC)"
	@echo -e "$(PRTP)CXX                \e[0m = $(CXX)"
	@echo -e "$(PRTP)FC                 \e[0m = $(FC)"
	@echo -e "$(PRTP)MAKEDEPEND         \e[0m = $(MAKEDEPEND)"
	@echo -e "$(PRTP)LD                 \e[0m = $(LD)"
	@echo -e "$(PRTP)SHARED_LD          \e[0m = $(SHARED_LD)"
	@echo -e "$(PRTP)DYNAMIC_LD         \e[0m = $(DYNAMIC_LD)"
	@echo -e "$(PRTP)FFLAGS             \e[0m = $(FFLAGS)"
	@echo -e "$(PRTP)CFLAGS             \e[0m = $(CFLAGS)"
	@echo -e "$(PRTP)CXXFLAGS           \e[0m = $(CXXFLAGS)"
	@echo -e "$(PRTP)LDFLAGS            \e[0m = $(LDFLAGS)"
	@echo -e "$(PRTP)SHARED_LD_FLAGS    \e[0m = $(SHARED_LD_FLAGS)"
	@echo -e "$(PRTP)DYNAMIC_LD_FLAGS   \e[0m = $(DYNAMIC_LD_FLAGS)"
	@echo -e "$(PRTP)LDPATH             \e[0m = $(LDPATH)"
	@echo -e "$(PRTP)RPATH              \e[0m = $(RPATH)"
	@echo -e "$(PRTP)DEFINES            \e[0m = $(DEFINES)"
	@echo -e "$(PRTP)FDEFINES           \e[0m = $(FDEFINES)"
	@echo -e "$(PRTP)LOCI_INCLUDES      \e[0m = $(LOCI_INCLUDES)"
	@echo -e "$(PRTP)INCLUDES           \e[0m = $(INCLUDES)"
	@echo -e "$(PRTP)LDLIBS             \e[0m = $(LDLIBS)"
	@echo -e "$(PRTP)LOCAL_LIBS         \e[0m = $(LOCAL_LIBS)"
	@echo -e "$(PRTP)LOCI_DEBUG         \e[0m = $(LOCI_DEBUG)"


