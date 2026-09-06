include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

DRIVER = ../FaceTransition/test_face_transition.exe
SOURCE_GRID = mixed.ugrid.txt
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

mixed.ugrid: $(SOURCE_GRID)
	cp $< $@

mixed.vog: mixed.ugrid
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m mixed \
	  > mixed.convert 2>&1

# Exercise the public transition handoff with adjacent hexahedral and prism
# roots. Both root topologies must retain valid face, cell, and node ancestry.
TestResults: $(DRIVER) mixed.vog FRC
	$(SERIAL) $(DRIVER) --scheduleoutput mixed > serial.log 2>&1 || { \
	  status=$$?; cat serial.log; exit $$status; }
	$(TWO_RANKS) $(DRIVER) --scheduleoutput mixed > parallel_2.log 2>&1 || { \
	  status=$$?; cat parallel_2.log; exit $$status; }
	$(THREE_RANKS) $(DRIVER) --scheduleoutput mixed > parallel_3.log 2>&1 || { \
	  status=$$?; cat parallel_3.log; exit $$status; }
	@echo "FVMAdapt2Tests/Module/MixedTransition: PASSED!" > $@
	@cat $@

clean:
	rm -rf output debug
	rm -f TestResults mixed.ugrid mixed.vog mixed.convert
	rm -f serial.log parallel_2.log parallel_3.log
	rm -f currentPlan nextPlan *.plan *.quality

distclean: clean

FRC:
