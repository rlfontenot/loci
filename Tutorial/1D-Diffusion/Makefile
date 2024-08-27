LOCI_BASE?=../..
LOCI_SRC?=../..

include $(LOCI_BASE)/Loci.conf


DIR_NAME	=	1D-Diffusion
TARGET1		=	heat
TARGET2		=	heat2

OBJ1		=	heat.o \
				main.o
OBJ2		=	heat2.o \
				main.o

################################################################################
# No changes required below this line
################################################################################
LOCI_FILES		=	$(wildcard *.loci)
LOCI_LPP_FILES	=	$(LOCI_FILES:.loci=.cc)
OBJS			=	$(OBJ1) $(OBJ2)
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))
JUNK			=	$(wildcard *~)  core make.depend make.depend.bak debug

default:	$(TARGET1)
all:		$(TARGET1) $(TARGET2)

$(TARGET1): $(OBJ1)
	@$(ECHO) "$(PRT) $(DIR_NAME)_$@ $(NC)"
	$(LD) -o $@ $(OBJ1) $(LOCAL_LIBS) $(LDFLAGS) $(LDLIBS)

$(TARGET2): $(OBJ2)
	@$(ECHO) "$(PRT) $(DIR_NAME)_$@ $(NC)"
	$(LD) -o $@ $(OBJ2) $(LOCAL_LIBS) $(LDFLAGS) $(LDLIBS)

clean:
	@$(ECHO) "$(PRT) $@_$(DIR_NAME) $(NC)"
	rm -rf  $(OBJS) $(JUNK) $(TARGET1) $(TARGET2)


distclean: clean
	@$(ECHO) "$(PRT) $@_$(DIR_NAME) $(NC)"
	rm -fr  $(LOCI_LPP_FILES) $(DEPEND_FILES)

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif