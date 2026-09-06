include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

INCLUDES += -I$(TEST_BASE)/../src/FVMAdapt2 \
            -I$(TEST_BASE)/../src/FVMAdapt2/library

TARGET = test_general_cell_paths.exe
FILES = test_general_cell_paths.cc
OBJS = $(call loci_file2objs,$(FILES))

SOURCE_GRID = ../General/tet.ugrid
LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH = $(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
SERIAL = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(SERIALRUN)
TWO_RANKS = LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) $(MPI_RUN) -np 2

.PHONY: default clean distclean FRC

default: TestResults

$(TARGET): $(OBJS) $(LOCI_BASE)/lib/libfvmadapt2func.$(LIB_SUFFIX)
	$(LD) -o $@ $(OBJS) -L$(LOCI_BASE)/lib -lfvmadapt2func $(LIBS) $(LDFLAGS)

tet.ugrid: $(SOURCE_GRID)
	cp $< $@

tet.vog: tet.ugrid
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m tet \
	  > tet.convert 2>&1

# Validate general-cell plans and query the scheduler-visible paths for an
# unrefined tetrahedron, a root split, and a nested split. MPI also covers an
# empty cell partition.
TestResults: $(TARGET) tet.vog FRC
	$(SERIAL) ./$(TARGET) tet > serial.log 2>&1 || { \
	  status=$$?; cat serial.log; exit $$status; }
	$(TWO_RANKS) ./$(TARGET) tet > parallel_2.log 2>&1 || { \
	  status=$$?; cat parallel_2.log; exit $$status; }
	@echo "FVMAdapt2Tests/Module/GeneralCellPaths: PASSED!" > $@
	@cat $@

clean:
	rm -rf output debug
	rm -f TestResults $(TARGET) $(OBJS) $(DEPEND_FILES)
	rm -f tet.ugrid tet.vog tet.convert serial.log parallel_2.log

distclean: clean

FRC:

DEPEND_FILES = $(OBJS:.o=.d)

ifeq ($(filter $(MAKECMDGOALS),clean distclean),)
-include $(DEPEND_FILES)
endif
