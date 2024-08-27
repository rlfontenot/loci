
LOCI_BASE?=..
LOCI_SRC?=..
include $(LOCI_SRC)/Loci.conf

TARGET	=	2dgv


#############################################################################
# No need to change rules below this line
#############################################################################
FILES			=	$(wildcard *.cpp)
OBJS			=	$(subst .cpp,.o,$(FILES))
DEPEND_FILES	=	$(subst .o,.d,$(OBJS))
JUNK			=	*~ core ii_files ti_files KCC_files so_locations rii_files
X11BASE			=	/usr/X11R6
CFLAGS			=	-O1
CXXFLAGS		=	-O1
LIBS			=	-lm -L$(X11BASE)/lib -lXt -lXaw -lX11
INCLUDES		=	-I$(LOCI_SRC)/include/2dgv \
					-I$(X11BASE)/include

# rule for generating dependencies from source files
%.d: %.c
	set -e; $(CC) -M $(CFLAGS) $(INCLUDES) $< \
	| sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' > $@; \
		[ -s $@ ] || rm -f $@
%.d: %.C
	set -e; $(CXX) -M $(CXXFLAGS) $(INCLUDES) $< \
	| sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' > $@; \
		[ -s $@ ] || rm -f $@
%.d: %.cc
	set -e; $(CXX) -M $(CXXFLAGS) $(INCLUDES) $< \
	| sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' > $@; \
		[ -s $@ ] || rm -f $@


default:	$(TARGET)
all:		$(TARGET)
spotless:	distclean clean_install

$(TARGET): $(OBJS)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(CXX) -o $@ $(OBJS) $(LIBS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

install: $(LOCI_BASE)/bin/. $(TARGET)
	@$(ECHO) "$(PRT) $@_$(TARGET) $(NC)"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	$(COPY) -u $(LOCI_SRC)/bin/$(TARGET) $(LOCI_BASE)/bin/$(TARGET)
endif

clean_install:
	@$(ECHO) "$(PRT) $@_$(TARGET) $(NC)"
	rm -rf $(LOCI_BASE)/bin/$(TARGET)

clean:
	@$(ECHO) "$(PRT) $@_$(TARGET) $(NC)"
	rm -fr $(LOCI_SRC)/bin/$(TARGET) $(TARGET) $(OBJS) $(JUNK)

distclean: clean
	@$(ECHO) "$(PRT) $@_$(TARGET) $(NC)"
	rm -fr $(DEPEND_FILES)

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif


