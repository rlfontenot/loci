LOCI_BASE?=../..
LOCI_SRC?=../..

include $(LOCI_SRC)/Loci.conf

MODULE_NAME	=	Datatypes
TARGETS		=	entities   \
				containers \
				relations  \
				helpers

################################################################################
# No changes required below this line
################################################################################
FILES			=	$(wildcard *.cpp)
OBJS			=	$(subst .cpp,.o,$(FILES))
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))

default:	$(TARGETS)
all:		$(TARGETS)

%: %.o
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $^ $(LDPATH) $(LDFLAGS) $(LDLIBS)

clean:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(OBJS) $(TARGETS)  $(JUNK)

distclean: clean
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(DEPEND_FILES)


ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif
