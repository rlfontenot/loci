###############################################################################
#
# Copyright 2008, Mississippi State University
#
# This file is part of the Loci Framework.
#
# The Loci Framework is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The Loci Framework is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
#
###############################################################################
LOCI_BASE?=..
LOCI_SRC?=..

include $(LOCI_SRC)/Loci.conf

MODULE_NAME	=	FVMtools

TARGET1		=	extract
TARGET2		=	plot3d2vog
TARGET3		=	make_periodic
TARGET4		=	ugrid2vog
TARGET5		=	cfd++2vog
TARGET6		=	fluent2vog
TARGET7		=	cobalt2vog
TARGET8		=	vogmerge
TARGET9		=	vogcut
TARGET10	=	vogcheck
TARGET11	=	extruder
TARGET12	=	vog2surf
TARGET13	=	ccm2vog
TARGET14	=	refmesh
TARGET15	=	marker
TARGET16	=	refine
TARGET17	=	cgns2vog
TARGET18	=	ugrid2cgns
TARGET19	=	cgns2ugrid
TARGET20	=	cgns2ensight
TARGET21	=	cgns2surf
TARGET22	=	extract_chem
TARGET23	=	ugridmirror
TARGET24	=	msh2vog
TARGET25	=	ccm2vogp

OBJ1	:=	extract.o         \
			extract_2dgv.o    \
			extract_ascii.o   \
			extract_cgns.o    \
			extract_cut.o     \
			extract_ensight.o \
			extract_fv.o      \
			extract_mean.o    \
			extract_surf.o    \
			extract_tecplot.o \
			extract_vtk.o     \
			particlePart.o    \
			surfacePart.o     \
			volumePart.o

OBJ2	:=	plot3d2vog.o vogtools.o
OBJ3	:=	make_periodic.o predicates.o vogtools.o
OBJ4	:=	ugrid2vog.o vogtools.o
OBJ5	:=	cfd++2vog.o vogtools.o
OBJ6	:=	fluent2vog.o vogtools.o
OBJ7	:=	cobalt2vog.o vogtools.o
OBJ8	:=	vogmerge.o
OBJ9	:=	vogcut.o vogtools.o
OBJ10	:=	vogCheck.o
OBJ11	:=	extruder.o
OBJ12	:=	vog2surf.o vogtools.o
OBJ13	:=	ccm2vog.o vogtools.o
OBJ14	:=	refmesh.o
OBJ15	:=	marker.o
OBJ16	:=	refine.o
OBJ17	:=	cgns2vog.o hash.o vogtools.o
OBJ18	:=	ugrid2cgns.o hash.o binaryio.o getargs.o
OBJ19	:=	cgns2ugrid.o hash.o binaryio.o getargs.o
OBJ20	:=	cgns2ensight.o
OBJ21	:=	cgns2surf.o hash.o
OBJ22	:=	extract_chem.o
OBJ23	:=	ugridmirror.o vogtools.o
OBJ24	:=	msh2vog.o vogtools.o
OBJ25	:=	ccm2vogp.o vogtools.o

DEPS	=	$(LOCI_SRC)/lib/libLoci.$(LIB_SUFFIX) \
			$(LOCI_SRC)/lib/libTools.$(LIB_SUFFIX) \
			$(LOCI_SRC)/lib/libsprng.$(LIB_SUFFIX)

DEPS10	=	$(DEPS) $(LOCI_SRC)/bin/lpp
DEPS13	=	$(DEPS) $(LOCI_SRC)/lib/libadf.a
DEPS14	=	$(DEPS) $(LOCI_SRC)/lib/fvmadapt_m.$(LIB_SUFFIX) $(LOCI_SRC)/lib/libfvmadaptfunc.$(LIB_SUFFIX)
DEPS15	=	$(DEPS) $(LOCI_SRC)/lib/fvmadapt_m.$(LIB_SUFFIX) $(LOCI_SRC)/lib/libfvmadaptfunc.$(LIB_SUFFIX)
DEPS16	=	$(DEPS) $(LOCI_SRC)/lib/fvmadapt_m.$(LIB_SUFFIX) $(LOCI_SRC)/lib/libfvmadaptfunc.$(LIB_SUFFIX)
DEPS25	=	$(DEPS) $(LOCI_SRC)/lib/libadf.a

ALL_TARGETS	=	$(TARGET1)  \
				$(TARGET2)  \
				$(TARGET3)  \
				$(TARGET4)  \
				$(TARGET5)  \
				$(TARGET6)  \
				$(TARGET7)  \
				$(TARGET8)  \
				$(TARGET9)  \
				$(TARGET10) \
				$(TARGET11) \
				$(TARGET12) \
				$(TARGET13) \
				$(TARGET14) \
				$(TARGET15) \
				$(TARGET16) \
				$(TARGET17) \
				$(TARGET18) \
				$(TARGET19) \
				$(TARGET20) \
				$(TARGET21) \
				$(TARGET23) \
				$(TARGET24) \
				$(TARGET25)

ALL_OBJ	=	$(OBJ1)  \
			$(OBJ2)  \
			$(OBJ3)  \
			$(OBJ4)  \
			$(OBJ5)  \
			$(OBJ6)  \
			$(OBJ7)  \
			$(OBJ8)  \
			$(OBJ9)  \
			$(OBJ10) \
			$(OBJ11) \
			$(OBJ12) \
			$(OBJ13) \
			$(OBJ14) \
			$(OBJ15) \
			$(OBJ16) \
			$(OBJ17) \
			$(OBJ18) \
			$(OBJ19) \
			$(OBJ20) \
			$(OBJ21) \
			$(OBJ23) \
			$(OBJ24) \
			$(OBJ25)

################################################################################
# No changes required below this line
################################################################################
LOCI_FILES		=	$(wildcard *.loci)
LOCI_LPP_FILES	=	$(LOCI_FILES:.loci=.cc)
JUNK			=	$(wildcard *~)  core
DEPEND_FILES	=	$(subst .o,.d,$(ALL_OBJ))
CFLAGS			+=	-O1
CXXFLAGS		+=	-O1
INCLUDES		=	-I$(LOCI_SRC)/include           \
					-I$(LOCI_SRC)/include/FVMtools  \
					-I$(LOCI_SRC)/include/libadf

default:	$(ALL_TARGETS)
spotless:	distclean clean_install

$(LOCI_LPP_FILES): $(LOCI_SRC)/bin/lpp

$(LOCI_SRC)/bin/lpp:
	$(MAKE) -C $(LOCI_SRC)/lpp
$(LOCI_SRC)/lib/libLoci.$(LIB_SUFFIX):
	$(MAKE) -C $(LOCI_SRC)/System
$(LOCI_SRC)/lib/libTools.$(LIB_SUFFIX):
	$(MAKE) -C $(LOCI_SRC)/Tools
$(LOCI_SRC)/lib/libsprng.$(LIB_SUFFIX):
	$(MAKE) -C $(LOCI_SRC)/sprng
$(LOCI_SRC)/lib/fvmadapt_m.$(LIB_SUFFIX):
	$(MAKE) -C $(LOCI_SRC)/FVMAdapt
$(LOCI_SRC)/lib/libfvmadaptfunc.$(LIB_SUFFIX):
	$(MAKE) -C $(LOCI_SRC)/FVMAdapt
$(LOCI_SRC)/lib/libadf.a:
	$(MAKE) -C libadf all

$(TARGET1): $(DEPS) $(OBJ1)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ1) $(LOCAL_LIBS) $(LDFLAGS) $(LDLIBS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET2): $(DEPS) $(OBJ2)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ2) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET3): $(DEPS) $(OBJ3)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ3) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET4): $(DEPS) $(OBJ4)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ4) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET5): $(DEPS) $(OBJ5)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ5) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET6): $(DEPS) $(OBJ6)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ6) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET7): $(DEPS) $(OBJ7)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ7) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET8): $(DEPS) $(OBJ8)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ8) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET9): $(DEPS) $(OBJ9)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ9) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET10): $(DEPS) $(OBJ10)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ10) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET11): $(DEPS) $(OBJ11)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ11) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET12): $(DEPS) $(OBJ12)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ12) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET13): LOCAL_LIBS+=libadf/libadf.a
$(TARGET13): $(DEPS13) $(OBJ13)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ13) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET14): LOCAL_LIBS+=-lfvmadaptfunc
$(TARGET14): $(DEPS14) $(OBJ14)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@  $(OBJ14) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET15): LOCAL_LIBS+=-lfvmadaptfunc
$(TARGET15): $(DEPS15) $(OBJ15)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@  $(OBJ15) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET16): LOCAL_LIBS+=-lfvmadaptfunc
$(TARGET16): $(DEPS16) $(OBJ16)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@  $(OBJ16) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET17): $(DEPS) $(OBJ17)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ17) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET18): $(DEPS) $(OBJ18)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ18) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET19): $(DEPS) $(OBJ19)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ19) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET20): $(DEPS) $(OBJ20)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ20)  $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET21): $(DEPS) $(OBJ21)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ21) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET22): $(DEPS) $(OBJ22)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ22) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET23): $(DEPS) $(OBJ23)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ23) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET24): $(DEPS) $(OBJ24)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ24) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

$(TARGET25): LOCAL_LIBS+=libadf/libadf.a
$(TARGET25): $(DEPS25) $(OBJ25)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(LD) -o $@ $(OBJ25) $(LOCAL_LIBS) $(LDLIBS) $(LDFLAGS)
	$(COPY) -u $@ $(LOCI_SRC)/bin/$@

install: $(ALL_TARGETS)
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	$(MAKE) -C libadf $@
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	for d in $(ALL_TARGETS); do $(COPY) $(LOCI_SRC)/bin/$$d $(LOCI_BASE)/bin/$$d; done
	$(COPY) extract_movie $(LOCI_BASE)/bin/extract_movie
endif

clean_install:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	for d in $(ALL_TARGETS); do rm -rf $(LOCI_BASE)/bin/$$d; done
	rm -rf $(LOCI_BASE)/bin/extract_movie
	$(MAKE) -C libadf $@

clean:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	@for d in $(ALL_TARGETS); do rm -rf $$d $(LOCI_SRC)/bin/$$d; done
	rm -rf  $(ALL_OBJ) $(JUNK)
	$(MAKE) -C libadf $@

distclean: clean
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(DEPEND_FILES) $(LOCI_LPP_FILES)
	$(MAKE) -C libadf $@

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif
