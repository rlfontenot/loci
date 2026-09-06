###############################################################################
#
# Copyright 2008-2026, Mississippi State University
#
# This file is part of the Loci Framework.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
###############################################################################

TESTS = Adaptation Transfer RemapPlan NodeTransition CellRemap
TEST_RESULTS = $(addsuffix /TestResults,$(TESTS))

.PHONY: default $(TESTS) clean distclean FRC

default: TestResults

TestResults: $(TESTS)
	cat $(TEST_RESULTS) > $@
	rm $(TEST_RESULTS)
	@cat $@

$(TESTS): FRC
	$(MAKE) -C $@ LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

clean:
	rm -f TestResults
	@for test in $(TESTS); do \
	  $(MAKE) -C $$test LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean || exit $$?; \
	done

distclean:
	rm -f TestResults
	@for test in $(TESTS); do \
	  $(MAKE) -C $$test LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean || exit $$?; \
	done

FRC:
