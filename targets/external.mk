################################################################################
# Author:      Mark A. Hunt (CFDRC)
# Date:        2024-07-04
# Description: This is the "external" targets filefile that holds the loci models.
################################################################################
include $(LOCI_SRC)/sys.conf

ifeq ($(shell uname -s),Darwin)
COPY	=	rsync -L
ESC		=	\\033[
ECHO	=	echo
SED_I	=	sed -i ''
else
COPY	=	cp
ESC		=	\e[
ECHO	=	echo -e
SED_I	=	sed -i
endif
PRT		=	$(ESC)1;33;44m
NC		=	$(ESC)0m

GKlib:
	@$(ECHO) "$(PRT) $@ $(NC)"
ifeq ($(INSTALL_GKLIB),1)
	$(LOCI_SRC)/targets/deps/install_$@.sh $(LOCI_SRC) $(GKLIB_BASE)
endif


METIS: GKlib
	@$(ECHO) "$(PRT) $@ $(NC)"
ifeq ($(INSTALL_METIS),1)
	$(LOCI_SRC)/targets/deps/install_$@.sh $(LOCI_SRC) $(METIS_BASE) $(GKLIB_BASE)
endif


ParMETIS: METIS
	@$(ECHO) "$(PRT) $@ $(NC)"
ifeq ($(INSTALL_PARMETIS),1)
	$(LOCI_SRC)/targets/deps/install_$@.sh $(LOCI_SRC) $(PARMETIS_BASE) $(GKLIB_BASE) \
		$(METIS_BASE)
endif


clean_GKlib:
	@$(ECHO) "$(PRT) $@ $(NC)"
	rm -rf $(LOCI_SRC)/$(shell echo $@ | cut -d"_" -f2)
clean_METIS:
	@$(ECHO) "$(PRT) $@ $(NC)"
	rm -rf $(LOCI_SRC)/$(shell echo $@ | cut -d"_" -f2)
clean_ParMETIS:
	@$(ECHO) "$(PRT) $@ $(NC)"
	rm -rf $(LOCI_SRC)/$(shell echo $@ | cut -d"_" -f2)

distclean_GKlib:
	$(MAKE) clean_$(shell echo $@ | cut -d"_" -f2)
distclean_METIS:
	$(MAKE) clean_$(shell echo $@ | cut -d"_" -f2)
distclean_ParMETIS:
	$(MAKE) clean_$(shell echo $@ | cut -d"_" -f2)

install_GKlib: GKlib
	@$(ECHO) "$(PRT) $@ $(NC)"
ifeq ($(INSTALL_GKLIB),1)
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	mkdir -p -- "$(LOCI_BASE)"
	$(COPY) -r $(LOCI_SRC)/$(shell echo $@ | cut -d"_" -f2) \
	           $(LOCI_BASE)/
endif
endif

install_METIS: METIS
	@$(ECHO) "$(PRT) $@ $(NC)"
ifeq ($(INSTALL_METIS),1)
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	mkdir -p -- "$(LOCI_BASE)"
	$(COPY) -r $(LOCI_SRC)/$(shell echo $@ | cut -d"_" -f2) \
	           $(LOCI_BASE)/
endif
endif

install_ParMETIS: ParMETIS
	@$(ECHO) "$(PRT) $@ $(NC)"
ifeq ($(INSTALL_PARMETIS),1)
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	mkdir -p -- "$(LOCI_BASE)"
	$(COPY) -r $(LOCI_SRC)/$(shell echo $@ | cut -d"_" -f2) \
	           $(LOCI_BASE)/
endif
endif

clean_install_GKlib:
	@$(ECHO) "$(PRT) $@ $(NC)"
ifeq ($(INSTALL_GKLIB),1)
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/$(shell echo $@ | cut -d"_" -f3)
endif
endif

clean_install_METIS:
	@$(ECHO) "$(PRT) $@ $(NC)"
ifeq ($(INSTALL_METIS),1)
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/$(shell echo $@ | cut -d"_" -f3)
endif
endif

clean_install_ParMETIS:
	@$(ECHO) "$(PRT) $@ $(NC)"
ifeq ($(INSTALL_PARMETIS),1)
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/$(shell echo $@ | cut -d"_" -f3)
endif
endif


spotless_GKlib:
	$(MAKE) distclean_$(shell     echo $@ | cut -d"_" -f2)
	$(MAKE) clean_install_$(shell echo $@ | cut -d"_" -f2)
spotless_METIS:
	$(MAKE) distclean_$(shell     echo $@ | cut -d"_" -f2)
	$(MAKE) clean_install_$(shell echo $@ | cut -d"_" -f2)
spotless_ParMETIS:
	$(MAKE) distclean_$(shell     echo $@ | cut -d"_" -f2)
	$(MAKE) clean_install_$(shell echo $@ | cut -d"_" -f2)
