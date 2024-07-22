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
	@echo -e "\e[$(PRT) $(DIR_NAME)_$@ \e[0m"
	$(LD) -o $@ $(OBJ1) $(LOCAL_LIBS) $(LDFLAGS) $(LDLIBS)

$(TARGET2): $(OBJ2)
	@echo -e "\e[$(PRT) $(DIR_NAME)_$@ \e[0m"
	$(LD) -o $@ $(OBJ2) $(LOCAL_LIBS) $(LDFLAGS) $(LDLIBS)

clean:
	@echo -e "\e[$(PRT) $@_$(DIR_NAME) \e[0m"
	rm -rf  $(OBJS) $(JUNK) $(TARGET1) $(TARGET2)


distclean: clean
	@echo -e "\e[$(PRT) $@_$(DIR_NAME) \e[0m"
	rm -fr  $(LOCI_LPP_FILES) $(DEPEND_FILES)

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif