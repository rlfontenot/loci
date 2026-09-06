include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

INCLUDES += -I$(TEST_BASE)/../src/FVMAdapt2/library \
            -I$(TEST_BASE)/contrib/doctest

TARGET = test_face_transition.exe
FILES = test_face_transition.loci
OBJS = $(call loci_file2objs,$(FILES))
LOCI_INTERMEDIATE_FILES = $(FILES:.loci=.$(LPP_I_SUFFIX))

SOURCE_GRID = ../Hex/hex.ugrid
LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
SERIAL = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(SERIALRUN)
TWO_RANKS = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(MPI_RUN) -np 2
THREE_RANKS = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(MPI_RUN) -np 3

.PHONY: default clean distclean FRC

default: TestResults

$(TARGET): $(OBJS) $(LOCI_BASE)/lib/libfvmadapt2func.$(LIB_SUFFIX)
	$(LD) -o $@ $(OBJS) -L$(LOCI_BASE)/lib -lfvmadapt2func $(LIBS) $(LDFLAGS)

hex.ugrid: $(SOURCE_GRID)
	cp $< $@

hex.vog: hex.ugrid
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m hex \
	  > hex.convert 2>&1

# Refine one hexahedral root, refine one existing leaf, and derefine while
# deliberately rebuilding the optional face-state cache from currentPlan.
# MPI covers repartitioned coarsening and empty local face partitions; serial
# also collapses to the root and verifies a different split direction receives
# distinct persistent child-cell and intra-cell face identities.
TestResults: $(TARGET) hex.vog FRC
	$(SERIAL) ./$(TARGET) --scheduleoutput hex > serial.log 2>&1 || { \
	  status=$$?; cat serial.log; exit $$status; }
	$(TWO_RANKS) ./$(TARGET) --scheduleoutput hex > parallel_2.log 2>&1 || { \
	  status=$$?; cat parallel_2.log; exit $$status; }
	$(THREE_RANKS) ./$(TARGET) --scheduleoutput hex > parallel_3.log 2>&1 || { \
	  status=$$?; cat parallel_3.log; exit $$status; }
	@echo "FVMAdapt2Tests/Module/FaceTransition: PASSED!" > $@
	@cat $@

clean:
	rm -rf output debug
	rm -f TestResults $(TARGET) $(OBJS) $(DEPEND_FILES) \
	  $(LOCI_INTERMEDIATE_FILES)
	rm -f hex.ugrid hex.vog hex.convert serial.log parallel_2.log parallel_3.log
	rm -f currentPlan nextPlan *.plan *.quality

distclean: clean
	rm -f $(OBJS:.o=.d)

FRC:

DEPEND_FILES = $(OBJS:.o=.d)

ifeq ($(filter $(MAKECMDGOALS),clean distclean),)
-include $(DEPEND_FILES)
endif
