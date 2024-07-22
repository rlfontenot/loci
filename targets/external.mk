include $(LOCI_SRC)/sys.conf
PRT	=	1;33;44m

GKlib:
	@echo -e "\e[$(PRT) $@ \e[0m"
ifeq ($(INSTALL_GKLIB),1)
	$(LOCI_SRC)/targets/deps/install_$@.sh $(LOCI_SRC) $(GKLIB_BASE)
endif


METIS: GKlib
	@echo -e "\e[$(PRT) $@ \e[0m"
ifeq ($(INSTALL_METIS),1)
	$(LOCI_SRC)/targets/deps/install_$@.sh $(LOCI_SRC) $(METIS_BASE) $(GKLIB_BASE)
endif


ParMETIS: METIS
	@echo -e "\e[$(PRT) $@ \e[0m"
ifeq ($(INSTALL_PARMETIS),1)
	$(LOCI_SRC)/targets/deps/install_$@.sh $(LOCI_SRC) $(PARMETIS_BASE) $(GKLIB_BASE) \
		$(METIS_BASE)
endif


clean_%:
	@echo -e "\e[$(PRT) $@ \e[0m"
	rm -rf $(LOCI_SRC)/$(*F)

distclean_%:
	$(MAKE) clean_$(*F)

install_%:
	@echo -e "\e[$(PRT) $@ \e[0m"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	cp -u $(LOCI_SRC)/$(*F) $(LOCI_BASE)/$(*F)
endif

clean_install_%:
	@echo -e "\e[$(PRT) $@ \e[0m"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/$(*F)
endif

clean_install_%:
	@echo -e "\e[$(PRT) $@ \e[0m"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/$(*F)
endif

spotless_%:
	$(MAKE) distclean_$(*F)
	$(MAKE) clean_install_$(*F)
