include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

DRIVER = ../FaceTransition/test_face_transition.exe
SOURCE_GRID = ../Prism/prism.ugrid
LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
SERIAL = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(SERIALRUN)
TWO_RANKS = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(MPI_RUN) -np 2
THREE_RANKS = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(MPI_RUN) -np 3

.PHONY: default clean distclean FRC

default: TestResults

$(DRIVER):
	$(MAKE) -C ../FaceTransition LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) \
	  test_face_transition.exe

prism.ugrid: $(SOURCE_GRID)
	cp $< $@

prism.vog: prism.ugrid
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m prism \
	  > prism.convert 2>&1

# Exercise the same pure-prism plan through initial and deeper refinement,
# partial derefinement, and cache-independent currentPlan replay. Boundary
# faces must retain boundary ancestry; created and removed faces are internal.
TestResults: $(DRIVER) prism.vog FRC
	$(SERIAL) $(DRIVER) --scheduleoutput prism > serial.log 2>&1 || { \
	  status=$$?; cat serial.log; exit $$status; }
	$(TWO_RANKS) $(DRIVER) --scheduleoutput prism > parallel_2.log 2>&1 || { \
	  status=$$?; cat parallel_2.log; exit $$status; }
	$(THREE_RANKS) $(DRIVER) --scheduleoutput prism > parallel_3.log 2>&1 || { \
	  status=$$?; cat parallel_3.log; exit $$status; }
	@echo "FVMAdapt2Tests/Module/PrismTransition: PASSED!" > $@
	@cat $@

clean:
	rm -rf output debug
	rm -f TestResults prism.ugrid prism.vog prism.convert
	rm -f serial.log parallel_2.log parallel_3.log
	rm -f currentPlan nextPlan *.plan *.quality

distclean: clean

FRC:
