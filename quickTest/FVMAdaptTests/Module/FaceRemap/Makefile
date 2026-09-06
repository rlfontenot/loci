include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf

INCLUDES += -I$(TEST_BASE)/contrib/doctest
TARGET = test_face_remap.exe
LIB_PATH = $(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
RUN = LD_LIBRARY_PATH=$(LIB_PATH) $(MPI_RUN)

.PHONY: default clean distclean FRC
default: TestResults

$(TARGET): test_face_remap.o $(LOCI_BASE)/lib/libfvmadapt2func.$(LIB_SUFFIX)
	$(LD) -o $@ test_face_remap.o -L$(LOCI_BASE)/lib -lfvmadapt2func $(LIBS) $(LDFLAGS)

tet.vog: tet.ugrid
	$(RUN) -np 1 $(LOCI_BASE)/bin/ugrid2vog -o -m tet > tet.convert 2>&1

TestResults: $(TARGET) tet.vog FRC
	@for ranks in 1 2 3; do \
	  $(RUN) -np $$ranks ./$(TARGET) --scheduleoutput > mpi_$$ranks.log 2>&1 || { \
	    status=$$?; cat mpi_$$ranks.log; exit $$status; }; \
	done
	@echo "FVMAdaptTests/Module/FaceRemap: PASSED!" > $@
	@cat $@

clean:
	rm -f $(TARGET) test_face_remap.o test_face_remap.d TestResults
	rm -f tet.vog tet.convert mpi_*.log currentPlan nextPlan *.plan *.quality
	rm -rf output debug

distclean: clean
FRC:

ifeq ($(filter $(MAKECMDGOALS),clean distclean),)
-include test_face_remap.d
endif
