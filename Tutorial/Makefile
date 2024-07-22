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

MODULE_NAME	=	docs
DIRS		=	1D-Diffusion heat Datatypes lfvm life

default:
	@echo -e "\e[$(PRT) $(MODULE_NAME) \e[0m"
	$(MAKE) -C docs

spotless:	distclean clean_install


install: default $(LOCI_BASE)/docs/.
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	for d in $(DIRS); do mkdir -p $(LOCI_BASE)/docs/$$d; done
	cp -u $(LOCI_SRC)/Tutorial/docs/tutorial.pdf        $(LOCI_BASE)/docs/
	cp -u $(LOCI_SRC)/Tutorial/LociTutorialSlides/*.pdf $(LOCI_BASE)/docs/

	cp -u $(LOCI_SRC)/Tutorial/1D-Diffusion/Makefile $(LOCI_BASE)/docs/1D-Diffusion/
	cp -u $(LOCI_SRC)/Tutorial/1D-Diffusion/*.loci   $(LOCI_BASE)/docs/1D-Diffusion/
	cp -u $(LOCI_SRC)/Tutorial/1D-Diffusion/*.lh     $(LOCI_BASE)/docs/1D-Diffusion/
	cp -u $(LOCI_SRC)/Tutorial/1D-Diffusion/*.cpp    $(LOCI_BASE)/docs/1D-Diffusion/
	cp -u $(LOCI_SRC)/Tutorial/1D-Diffusion/*.vars   $(LOCI_BASE)/docs/1D-Diffusion/

	cp -u $(LOCI_SRC)/Tutorial/heat/Makefile         $(LOCI_BASE)/docs/heat/
	cp -u $(LOCI_SRC)/Tutorial/heat/*.loci           $(LOCI_BASE)/docs/heat/
	cp -u $(LOCI_SRC)/Tutorial/heat/*.lh             $(LOCI_BASE)/docs/heat/
	cp -u $(LOCI_SRC)/Tutorial/heat/*.vars           $(LOCI_BASE)/docs/heat/

	cp -u $(LOCI_SRC)/Tutorial/Datatypes/Makefile    $(LOCI_BASE)/docs/Datatypes/
	cp -u $(LOCI_SRC)/Tutorial/Datatypes/*.cpp       $(LOCI_BASE)/docs/Datatypes/

	cp -u $(LOCI_SRC)/Tutorial/lfvm/Makefile         $(LOCI_BASE)/docs/lfvm/
	cp -u $(LOCI_SRC)/Tutorial/lfvm/*.h              $(LOCI_BASE)/docs/lfvm/
	cp -u $(LOCI_SRC)/Tutorial/lfvm/*.cpp            $(LOCI_BASE)/docs/lfvm/

	cp -u $(LOCI_SRC)/Tutorial/life/Makefile         $(LOCI_BASE)/docs/life/
	cp -u $(LOCI_SRC)/Tutorial/life/*.cpp            $(LOCI_BASE)/docs/life/
endif

clean_install:
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"
	for d in $(DIRS); do rm -rf $(LOCI_BASE)/docs/$$d; done
	rm -rf $(LOCI_BASE)/docs/*.pdf

clean:
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"
	$(MAKE) -C docs $@

distclean: clean
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"



