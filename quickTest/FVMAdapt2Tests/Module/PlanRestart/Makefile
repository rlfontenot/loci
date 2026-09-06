include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

INCLUDES += -I$(TEST_BASE)/../src/FVMAdapt2/library

TARGET = test_plan_restart.exe
FILES = test_plan_restart.cc
OBJS = $(call loci_file2objs,$(FILES))

SERIAL_GRID = ../Hex/hex.ugrid
GENERAL_GRID = ../General/tet.ugrid
MPI_GRID = ../../../FVMModUnitTests/tests/block.vog
LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
SERIAL = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(SERIALRUN)
TWO_RANKS = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(MPI_RUN) -np 2
THREE_RANKS = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(MPI_RUN) -np 3

.PHONY: default clean distclean FRC

default: TestResults

$(TARGET): $(OBJS) $(LOCI_BASE)/lib/libfvmadapt2func.$(LIB_SUFFIX)
	$(LD) -o $@ $(OBJS) -L$(LOCI_BASE)/lib -lfvmadapt2func $(LIBS) $(LDFLAGS)

hex.ugrid: $(SERIAL_GRID)
	cp $< $@

hex.vog: hex.ugrid
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m hex \
	  > hex.convert 2>&1

block.vog: $(MPI_GRID)
	cp $< $@

tet.ugrid: $(GENERAL_GRID)
	cp $< $@

tet.vog: tet.ugrid
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m tet \
	  > tet.convert 2>&1

# Reconstruct a nontrivial accepted plan, install its canonical face and cell
# identities after redistribution, then use its compact face state as the
# exact source of a subsequent transition. MPI uses a 64-cell pure-hex grid so
# this test remains independent of the separate empty-partition I/O fix. The
# tetrahedron checks that face identities work without general-node ancestry.
TestResults: $(TARGET) hex.vog block.vog tet.vog FRC
	$(SERIAL) ./$(TARGET) hex > serial.log 2>&1 || { \
	  status=$$?; cat serial.log; exit $$status; }
	$(SERIAL) ./$(TARGET) --general tet \
	  > general.log 2>&1 || { \
	  status=$$?; cat general.log; exit $$status; }
	$(TWO_RANKS) ./$(TARGET) block > parallel_2.log 2>&1 || { \
	  status=$$?; cat parallel_2.log; exit $$status; }
	$(THREE_RANKS) ./$(TARGET) block > parallel_3.log 2>&1 || { \
	  status=$$?; cat parallel_3.log; exit $$status; }
	@echo "FVMAdapt2Tests/Module/PlanRestart: PASSED!" > $@
	@cat $@

clean:
	rm -rf output debug
	rm -f TestResults $(TARGET) $(OBJS) $(DEPEND_FILES)
	rm -f hex.ugrid hex.vog block.vog tet.ugrid tet.vog hex.convert \
	  tet.convert serial.log general.log parallel_2.log parallel_3.log
	rm -f plan_restart.h5 missing_plan_restart_weights.h5
	rm -f currentPlan nextPlan *.plan *.quality

distclean: clean

FRC:

DEPEND_FILES = $(OBJS:.o=.d)

ifeq ($(filter $(MAKECMDGOALS),clean distclean),)
-include $(DEPEND_FILES)
endif
