LOCI_BASE?=..
LOCI_SRC?=..

include $(LOCI_SRC)/Loci.conf
DIR_NAME	=	sprng
MODULE_NAME	=	$(DIR_NAME)
LIB_NAME1	=	lib$(MODULE_NAME).a
LIB_NAME2	=	lib$(MODULE_NAME).$(LIB_SUFFIX)

################################################################################
# No changes required below this line
################################################################################
FILES			=	$(wildcard *.c)
OBJS			=	$(subst .c,.o,$(FILES))
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))
JUNK			=	*~ core ii_files ti_files KCC_files so_locations rii_files
INCLUDES		=	-I$(LOCI_SRC)/include/sprng \
					-I$(LOCI_SRC)/include
DEFINES			+=	-DSPRNG_MPI -DAdd_
CFLAGS			+=	-O1

default:	$(LIB_NAME2)
all:		$(LIB_NAME1) $(LIB_NAME2)
spotless:	distclean clean_install

$(LIB_NAME1): $(OBJS)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(RANLIB) $@
	$(COPY) -u $@ $(LOCI_SRC)/lib/$@

$(LIB_NAME2): $(OBJS)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(DYNAMIC_LD) $(DYNAMIC_LD_FLAGS) $@ $(LDFLAGS) $(OBJS)
	$(COPY) -u $@ $(LOCI_SRC)/lib/$@

install: $(LOCI_BASE)/lib/. $(LIB_NAME2)
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	$(COPY) -u $(LOCI_SRC)/lib/$(LIB_NAME2) $(LOCI_BASE)/lib/$(LIB_NAME2)
endif

clean_install:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME1)
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME2)

clean:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(LOCI_SRC)/lib/$(LIB_NAME1) $(LOCI_SRC)/lib/$(LIB_NAME2)
	rm -fr $(LIB_NAME1) $(LIB_NAME2) $(OBJS) $(JUNK)

distclean: clean
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(DEPEND_FILES)

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif

