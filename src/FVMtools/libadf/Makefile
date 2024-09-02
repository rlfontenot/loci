LOCI_BASE?=../..
LOCI_SRC?=../..

include $(LOCI_SRC)/Loci.conf

DIR_NAME	=	libadf
MODULE_NAME	=	adf
LIB_NAME1	=	lib$(MODULE_NAME).a
LIB_NAME2	=	lib$(MODULE_NAME).$(LIB_SUFFIX)

################################################################################
# No changes required below this line
################################################################################
FILES			=	$(wildcard *.c)
OBJS			=	$(subst .c,.o,$(FILES))
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))
JUNK			=	*~ core ii_files ti_files so_locations rii_files
INCLUDES		=	-I$(LOCI_SRC)/include/libadf
CFLAGS			+=	-O1

default:	$(LIB_NAME1)
all:		$(LIB_NAME1) $(LIB_NAME2)
spotless:	distclean clean_install

$(LIB_NAME1): $(OBJS)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(RANLIB) $@
	$(COPY) $@ $(LOCI_SRC)/lib/$@

$(LIB_NAME2): $(OBJS)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(SHARED_LD) $(SHARED_LD_FLAGS) $@ $(LDFLAGS) $(OBJS)
	$(COPY) $@ $(LOCI_SRC)/lib/$@

install: $(LOCI_BASE)/lib/. $(LIB_NAME1) $(LIB_NAME2)
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	$(COPY) $(LOCI_SRC)/lib/$(LIB_NAME1) $(LOCI_BASE)/lib/$(LIB_NAME1)
	$(COPY) $(LOCI_SRC)/lib/$(LIB_NAME2) $(LOCI_BASE)/lib/$(LIB_NAME2)
endif

clean_install:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME1)
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME2)

clean:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(OBJS) $(LIB_NAME1) $(LIB_NAME2) $(JUNK)
	rm -fr $(LOCI_SRC)/lib/$(LIB_NAME1)
	rm -fr $(LOCI_SRC)/lib/$(LIB_NAME2)

distclean: clean
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(DEPEND_FILES)


#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif
