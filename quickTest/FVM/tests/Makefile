TEST_FILES=cylinder.test fixedMass.test  wedge.test shock_tube.test 


#include $(LOCI_BASE)/test.conf
include ../../test.conf
#include test.conf

default: TestResults

TestResults: $(TEST_FILES) FRC
	cat $(TEST_FILES) > TestResults; rm -f $(TEST_FILES)

FRC:

%.test : %.vars
	FLOWTEST="$(FLOWTEST)" EXTRACT="$(EXTRACT)" NDIFF=../ndiff ./tests.sh $*.vars

clean:
	rm -fr TEST_* TestResults $(TEST_FILES)

distclean:
	rm -fr TEST_* TestResults $(TEST_FILES)


