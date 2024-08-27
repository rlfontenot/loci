LOCI_BASE?=../..
LOCI_SRC?=../..

include $(LOCI_BASE)/Loci.conf


TARGET	=	heat

################################################################################
# No changes required below this line
################################################################################
LOCI_FILES		=	$(wildcard *.loci)
LOCI_LPP_FILES	=	$(LOCI_FILES:.loci=.cc)
OBJS			=	$(subst .loci,.o,$(LOCI_FILES))
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))
JUNK			=	$(wildcard *~)  core make.depend make.depend.bak debug output

default:	$(TARGET)
all:		$(TARGET)

$(TARGET): $(OBJS)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJS) $(LOCAL_LIBS) $(LDFLAGS) $(LDLIBS)

clean:
	@$(ECHO) "$(PRT) $@_$(TARGET) $(NC)"
	rm -rf  $(OBJS) $(JUNK) $(TARGET)


distclean: clean
	@$(ECHO) "$(PRT) $@_$(TARGET) $(NC)"
	rm -fr  $(LOCI_LPP_FILES) $(DEPEND_FILES)

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif