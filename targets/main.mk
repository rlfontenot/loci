################################################################################
# Author:      Mark A. Hunt (CFDRC)
# Date:        October 14, 2022
# Description: This is the "targets" file that holds the loci models.
################################################################################
LOCI_SRC	?=	$(shell pwd)/{{OBJDIR}}

################################################################################
# Wanted to make this the default when calling the `make` command from root
# directory.
################################################################################
all:
	$(MAKE) -C {{OBJDIR}} default

################################################################################
################################################################################
include {{OBJDIR}}/targets/external.mk

################################################################################
################### (1) These are the MAKE definitions #########################
################################################################################
.PHONY: 2dgv
2dgv:
	$(MAKE) -C {{OBJDIR}}/$@ $@

copy_config:
ifneq ($(LOCI_SRC)/{{OBJDIR}},$(LOCI_BASE))
	mkdir -p $(LOCI_BASE)
	$(COPY) {{OBJDIR}}/*.conf $(LOCI_BASE)/
	$(SED_I) "s,LOCI_SRC,LOCI_BASE,g" $(LOCI_BASE)/Loci.conf
	$(SED_I) "s,LOCI_SRC,LOCI_BASE,g" $(LOCI_BASE)/sys.conf
endif


devhelp:
	$(MAKE) -C {{OBJDIR}} $@

FVMAdapt:
	$(MAKE) -C {{OBJDIR}}/$@

.PHONY: FVMMod
FVMMod:
	$(MAKE) -C {{OBJDIR}}/$@

.PHONY: FVMOverset
FVMOverset:
	$(MAKE) -C {{OBJDIR}}/$@

.PHONY: FVMtools
FVMtools:
	$(MAKE) -C {{OBJDIR}}/$@

list:
	$(MAKE) -C {{OBJDIR}} $@

list_installed:
	$(MAKE) -C {{OBJDIR}} $@

.PHONY: lpp
lpp:
	$(MAKE) -C {{OBJDIR}}/$@

.PHONY: sprng
sprng:
	$(MAKE) -C {{OBJDIR}}/$@

.PHONY: System
System:
	$(MAKE) -C {{OBJDIR}}/$@

.PHONY: Tools
Tools:
	$(MAKE) -C {{OBJDIR}}/$@

.PHONY: Tutorial
Tutorial:
	$(MAKE) -C {{OBJDIR}}/$@

################################################################################
################### (2) These are the pattern recipes ##########################
################################################################################

install_%: copy_config
	$(MAKE) -C {{OBJDIR}}/$(*F) install
clean_install_%: copy_config
	$(MAKE) -C {{OBJDIR}}/$(*F) clean_install
clean_%:
	$(MAKE) -C {{OBJDIR}}/$(*F) clean
distclean_%:
	$(MAKE) -C {{OBJDIR}}/$(*F) distclean
spotless_%:
	$(MAKE) -C {{OBJDIR}}/$(*F) spotless

install_include: copy_config
	$(MAKE) -C {{OBJDIR}}/System $@

clean_install_include:
	$(MAKE) -C {{OBJDIR}}/System $@