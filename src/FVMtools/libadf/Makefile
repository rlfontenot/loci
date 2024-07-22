LOCI_BASE?=../..
LOCI_SRC?=../..

include $(LOCI_SRC)/Loci.conf

DIR_NAME	=	libadf
MODULE_NAME	=	adf
LIB_NAME1	=	lib$(MODULE_NAME).a
LIB_NAME2	=	lib$(MODULE_NAME).so

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
	@echo -e "\e[$(PRT) $@ \e[0m"
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(RANLIB) $@
	cp $@ $(LOCI_SRC)/lib/$@

$(LIB_NAME2): $(OBJS)
	@echo -e "\e[$(PRT) $@ \e[0m"
	$(SHARED_LD) $(SHARED_LD_FLAGS) $@ $(LDFLAGS) $(OBJS)
	cp $@ $(LOCI_SRC)/lib/$@

install: $(LOCI_BASE)/lib/. $(LIB_NAME1) $(LIB_NAME2)
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	cp $(LOCI_SRC)/lib/$(LIB_NAME1) $(LOCI_BASE)/lib/$(LIB_NAME1)
	cp $(LOCI_SRC)/lib/$(LIB_NAME2) $(LOCI_BASE)/lib/$(LIB_NAME2)
endif

clean_install:
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME1)
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME2)

clean:
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"
	rm -fr $(OBJS) $(LIB_NAME1) $(LIB_NAME2) $(JUNK)
	rm -fr $(LOCI_SRC)/lib/$(LIB_NAME1)
	rm -fr $(LOCI_SRC)/lib/$(LIB_NAME2)

distclean: clean
	@echo -e "\e[$(PRT) $@_$(MODULE_NAME) \e[0m"
	rm -fr $(DEPEND_FILES)


#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif
