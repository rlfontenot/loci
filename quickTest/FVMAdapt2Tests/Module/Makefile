TESTS = RefMesh UniformHex Prism Mixed PlanInputs ExtrudedModes \
	FaceRemap FaceTransition GeneralCellPaths MixedTransition PlanRestart \
	PrismTransition RefinementState
TEST_RESULTS = $(addsuffix /TestResults,$(TESTS))

.PHONY: default $(TESTS) clean distclean FRC

# Scenarios share input fixtures. Keep their execution ordered while allowing
# each scenario's build to use the caller's make job limit.
.NOTPARALLEL: TestResults

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
	@for test in $(TESTS); do \
	  $(MAKE) -C $$test LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean || exit $$?; \
	done

distclean:
	rm -f TestResults
	@for test in $(TESTS); do \
	  $(MAKE) -C $$test LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean || exit $$?; \
	done

FRC:
