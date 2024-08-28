LOCI_BASE?=../..
LOCI_SRC?=../..

include $(LOCI_BASE)/Loci.conf


TARGET	=	lfvm

################################################################################
# No changes required below this line
################################################################################
FILES			=	$(wildcard *.cpp)
OBJS			=	$(subst .cpp,.o,$(FILES))
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))
JUNK			=	$(wildcard *~)  core debug output
INCLUDES		+=	-I.

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
	rm -fr $(DEPEND_FILES)

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif