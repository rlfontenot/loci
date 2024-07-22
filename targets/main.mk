################################################################################
# Author:      Mark A. Hunt (CFDRC)
# Date:        October 14, 2022
# Version:     v0.1
# Description: This is the "targets" file that holds the loci models.
################################################################################
LOCI_SRC	?=	$(shell pwd)/{{OBJDIR}}

include {{OBJDIR}}/targets/external.mk

################################################################################
################### (1) These are the MAKE definitions #########################
################################################################################
copy_config:
ifneq ($(LOCI_SRC)/{{OBJDIR}},$(LOCI_BASE))
	mkdir -p $(LOCI_BASE)
	cp {{OBJDIR}}/*.conf $(LOCI_BASE)/
	sed -i "s,LOCI_SRC,LOCI_BASE,g" $(LOCI_BASE)/Loci.conf
	sed -i "s,LOCI_SRC,LOCI_BASE,g" $(LOCI_BASE)/sys.conf
endif

all:
	$(MAKE) -C {{OBJDIR}} default

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
################### (2) These are the INSTALL definitions ######################
################################################################################
install_FVMAdapt: copy_config
	$(MAKE) -C {{OBJDIR}}/FVMAdapt install
install_FVMMod: copy_config
	$(MAKE) -C {{OBJDIR}}/FVMMod install
install_FVMOverset: copy_config
	$(MAKE) -C {{OBJDIR}}/FVMOverset install
install_FVMtools: copy_config
	$(MAKE) -C {{OBJDIR}}/FVMtools install
install_include: copy_config
	$(MAKE) -C {{OBJDIR}}/System $@
install_lpp: copy_config
	$(MAKE) -C {{OBJDIR}}/lpp install
install_sprng: copy_config
	$(MAKE) -C {{OBJDIR}}/sprng install
install_System: copy_config
	$(MAKE) -C {{OBJDIR}}/System install
install_Tools: copy_config
	$(MAKE) -C {{OBJDIR}}/Tools install
install_Tutorial: copy_config
	$(MAKE) -C {{OBJDIR}}/Tutorial install

################################################################################
################### (3) These are the CLEAN definitions ########################
################################################################################
clean_FVMAdapt:
	$(MAKE) -C {{OBJDIR}}/FVMAdapt clean
clean_FVMMod:
	$(MAKE) -C {{OBJDIR}}/FVMMod clean
clean_FVMOverset:
	$(MAKE) -C {{OBJDIR}}/FVMOverset clean
clean_FVMtools:
	$(MAKE) -C {{OBJDIR}}/FVMtools clean
clean_lpp:
	$(MAKE) -C {{OBJDIR}}/lpp clean
clean_sprng:
	$(MAKE) -C {{OBJDIR}}/sprng clean
clean_System:
	$(MAKE) -C {{OBJDIR}}/System clean
clean_Tools:
	$(MAKE) -C {{OBJDIR}}/Tools clean
clean_Tutorial:
	$(MAKE) -C {{OBJDIR}}/Tutorial clean

################################################################################
################### (4) These are the CLEAN_INSTALL definitions ################
################################################################################
clean_install_FVMAdapt:
	$(MAKE) -C {{OBJDIR}}/FVMAdapt clean_install
clean_install_FVMMod:
	$(MAKE) -C {{OBJDIR}}/FVMMod clean_install
clean_install_FVMOverset:
	$(MAKE) -C {{OBJDIR}}/FVMOverset clean_install
clean_install_FVMtools:
	$(MAKE) -C {{OBJDIR}}/FVMtools clean_install
clean_install_include:
	$(MAKE) -C {{OBJDIR}}/System $@
clean_install_lpp:
	$(MAKE) -C {{OBJDIR}}/lpp clean_install
clean_install_sprng:
	$(MAKE) -C {{OBJDIR}}/sprng clean_install
clean_install_System:
	$(MAKE) -C {{OBJDIR}}/System clean_install
clean_install_Tools:
	$(MAKE) -C {{OBJDIR}}/Tools clean_install
clean_install_Tutorial:
	$(MAKE) -C {{OBJDIR}}/Tutorial clean_install

################################################################################
################### (5) These are the DISTCLEAN definitions ####################
################################################################################
distclean_FVMAdapt:
	$(MAKE) -C {{OBJDIR}}/FVMAdapt distclean
distclean_FVMMod:
	$(MAKE) -C {{OBJDIR}}/FVMMod distclean
distclean_FVMOverset:
	$(MAKE) -C {{OBJDIR}}/FVMOverset distclean
distclean_FVMtools:
	$(MAKE) -C {{OBJDIR}}/FVMtools distclean
distclean_lpp:
	$(MAKE) -C {{OBJDIR}}/lpp distclean
distclean_sprng:
	$(MAKE) -C {{OBJDIR}}/sprng distclean
distclean_System:
	$(MAKE) -C {{OBJDIR}}/System distclean
distclean_Tools:
	$(MAKE) -C {{OBJDIR}}/Tools distclean
distclean_Tutorial:
	$(MAKE) -C {{OBJDIR}}/Tutorial distclean

################################################################################
################### (6) These are the SPOTLESS definitions #####################
################################################################################
spotless_FVMAdapt:
	$(MAKE) -C {{OBJDIR}}/FVMAdapt spotless
spotless_FVMMod:
	$(MAKE) -C {{OBJDIR}}/FVMMod spotless
spotless_FVMOverset:
	$(MAKE) -C {{OBJDIR}}/FVMOverset spotless
spotless_FVMtools:
	$(MAKE) -C {{OBJDIR}}/FVMtools spotless
spotless_lpp:
	$(MAKE) -C {{OBJDIR}}/lpp spotless
spotless_sprng:
	$(MAKE) -C {{OBJDIR}}/sprng spotless
spotless_System:
	$(MAKE) -C {{OBJDIR}}/System spotless
spotless_Tools:
	$(MAKE) -C {{OBJDIR}}/Tools spotless
spotless_Tutorial:
	$(MAKE) -C {{OBJDIR}}/Tutorial spotless
