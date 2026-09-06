TESTS = FaceTransition GeneralCellPaths MixedTransition PlanRestart \
	PrismTransition RefinementState
TEST_RESULTS = $(addsuffix /TestResults,$(TESTS))

.PHONY: default $(TESTS) clean distclean FRC

default: TestResults

TestResults: $(TESTS)
	cat $(TEST_RESULTS) > $@
	rm $(TEST_RESULTS)
	@cat $@

$(TESTS): FRC
	$(MAKE) -C $@ LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

# The mixed and prism qualifications reuse the generic face-transition driver.
MixedTransition: FaceTransition
PrismTransition: FaceTransition

clean:
	rm -f TestResults
	$(MAKE) -C FaceTransition LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C GeneralCellPaths LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C MixedTransition LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C PlanRestart LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C PrismTransition LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C RefinementState LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean

distclean:
	rm -f TestResults
	$(MAKE) -C FaceTransition LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C GeneralCellPaths LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C MixedTransition LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C PlanRestart LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C PrismTransition LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C RefinementState LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean

FRC:
