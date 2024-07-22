
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
	@echo -e "\e[$(PRT) $@ \e[0m"
	$(CXX) -o $@ $(OBJS) $(LIBS)
	cp -u $@ $(LOCI_SRC)/bin/$@

install: $(LOCI_BASE)/bin/. $(TARGET)
	@echo -e "\e[$(PRT) $@_$(TARGET) \e[0m"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	cp -u $(LOCI_SRC)/bin/$(TARGET) $(LOCI_BASE)/bin/$(TARGET)
endif

clean_install:
	@echo -e "\e[$(PRT) $@_$(TARGET) \e[0m"
	rm -rf $(LOCI_BASE)/bin/$(TARGET)

clean:
	@echo -e "\e[$(PRT) $@_$(TARGET) \e[0m"
	rm -fr $(LOCI_SRC)/bin/$(TARGET) $(TARGET) $(OBJS) $(JUNK)

distclean: clean
	@echo -e "\e[$(PRT) $@_$(TARGET) \e[0m"
	rm -fr $(DEPEND_FILES)

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif


