##############################################################################
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
##############################################################################
LOCI_BASE?=..
LOCI_SRC?=..
include $(LOCI_SRC)/Loci.conf

DIR_NAME	=	FVMAdapt
MODULE_NAME	=	fvmadapt
LIB_NAME1	=	$(MODULE_NAME)_m.so
LIB_NAME2	=	lib$(MODULE_NAME)func.$(LIB_SUFFIX)


OBJS1 :=	balance_general_cell.o          \
			balance_general_edge.o          \
			balance_general_edge_derefine.o \
			balance_general_face.o          \
			balance_general_face_derefine.o \
			balance_hex_cell.o              \
			balance_mixed_face.o            \
			balance_mixed_face_derefine.o   \
			balance_prism_cell.o            \
			balance_prism_face.o            \
			balance_prism_face_derefine.o   \
			balance_quad_face.o             \
			balance_quad_face_derefine.o    \
			cellplan_io.o                   \
			check_folded_face.o             \
			create_hexcell_maps.o           \
			create_prism_maps.o             \
			get_fine_face.o                 \
			get_fine_grid.o                 \
			get_hex_fine_face.o             \
			get_hex_fine_grid.o             \
			get_mixed_face_fine_grid.o      \
			get_node_remap.o                \
			get_prism_fine_face.o           \
			get_prism_fine_grid.o           \
			make_general_cellplan.o         \
			make_hex_cellplan.o             \
			make_prism_cellplan.o           \
			merge_general_edge.o            \
			merge_general_face.o            \
			merge_hex_face.o                \
			merge_mixed_face.o              \
			merge_prism_face.o              \
			parse_tree.o                    \
			quadface_rule.o                 \
			read_tag.o                      \
			reset_general_nums.o            \
			reset_hex_nums.o                \
			reset_offset.o                  \
			reset_prism_nums.o              \
			set_general_nums.o              \
			set_hex_nums.o                  \
			set_offset.o                    \
			set_prism_nums.o                \
			update_general_edge.o           \
			update_general_face.o           \
			update_mixed_face.o             \
			update_prism_face.o             \
			update_quad_face.o              \
			write_vog_module.o

OBJS2 := 	build_general_cell.o \
			build_hexcell.o      \
			build_prismcell.o    \
			classify_cell.o      \
			color_matrix.o       \
			diamondcell.o        \
			extract_hex_edge.o   \
			extract_hex_face.o   \
			extract_prism_face.o \
			face.o               \
			get_c1_general.o     \
			get_c1_hex.o         \
			get_c1_prism.o       \
			globals.o            \
			hexcell.o            \
			node_edge.o          \
			prism.o              \
			quadface.o           \
			read_par.o           \
			transfer_fc.o        \
			transfer_plan.o      \
			write_vog.o          \

################################################################################
# No changes required below this line
################################################################################
LOCI_FILES		:=	$(wildcard *.loci)
LOCI_LPP_FILES	:=	$(LOCI_FILES:.loci=.cc)
DEPEND_FILES	=	$(subst .o,.d,$(OBJS1),$(OBJS2))
JUNK			=	*~ core ti_files ii_files rii_files
INCLUDES		=	-I$(LOCI_SRC)/include/FVMAdapt \
					-I$(LOCI_SRC)/include          \
					$(PARMETIS_INCLUDE)            \
					$(LIBXML2_INCLUDE)
LOCAL_LIBS		=	$(LIBXML2_LDFLAGS) \
					$(LIBXML2_LDLIBS)  \
					-L$(LOCI_SRC)/lib -lLoci
FILE_DEPS		=	$(LOCI_SRC)/bin/lpp        \
					$(PARMETIS_BASE)/lib/libparmetis.$(LIB_SUFFIX)
LIB_DEPS		=	$(LOCI_SRC)/lib/libLoci.$(LIB_SUFFIX)

default:	$(LIB_NAME1) $(LIB_NAME2)
all:		$(LIB_NAME1) $(LIB_NAME2)
spotless:	distclean clean_install

$(LOCI_LPP_FILES): $(FILE_DEPS)
write_vog.d:       $(FILE_DEPS)

$(LOCI_SRC)/bin/lpp:
	$(MAKE) -C $(LOCI_SRC)/lpp
$(LOCI_SRC)/lib/libLoci.$(LIB_SUFFIX):
	$(MAKE) -C $(LOCI_SRC)/System
$(PARMETIS_BASE)/lib/libparmetis.$(LIB_SUFFIX):
	$(MAKE) -C $(LOCI_SRC) ParMETIS

$(LIB_NAME1): $(LIB_DEPS) $(OBJS1)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(SHARED_LD) $(SHARED_LD_FLAGS) $@ $(OBJS1) $(LOCAL_LIBS)
	$(COPY) -u $@ $(LOCI_SRC)/lib/$@

$(LIB_NAME2): $(LIB_DEPS) $(OBJS2)
	@$(ECHO) "$(PRT) $@ $(NC)"
	$(DYNAMIC_LD) $(DYNAMIC_LD_FLAGS) $@ $(OBJS2)
	$(COPY) -u $@ $(LOCI_SRC)/lib/$@


install: $(LOCI_BASE)/include/$(DIR_NAME)/. $(LOCI_BASE)/lib/.
install: $(LIB_NAME1) $(LIB_NAME2)
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	$(COPY) -u $(LOCI_SRC)/lib/$(LIB_NAME1) $(LOCI_BASE)/lib/$(LIB_NAME1)
	$(COPY) -u $(LOCI_SRC)/lib/$(LIB_NAME2) $(LOCI_BASE)/lib/$(LIB_NAME2)
	$(COPY) -u $(LOCI_SRC)/include/$(DIR_NAME)/*.h $(LOCI_BASE)/include/$(DIR_NAME)/
endif

clean_install:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME1)
	rm -rf $(LOCI_BASE)/lib/$(LIB_NAME2)
ifneq ($(LOCI_SRC),$(LOCI_BASE))
	rm -rf $(LOCI_BASE)/include/$(DIR_NAME)
endif

clean:
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(OBJS1) $(OBJS2) $(LIB_NAME1) $(LIB_NAME2)
	rm -fr $(LOCI_SRC)/lib/$(LIB_NAME1) $(LOCI_SRC)/lib/$(LIB_NAME2) $(JUNK)

distclean: clean
	@$(ECHO) "$(PRT) $@_$(MODULE_NAME) $(NC)"
	rm -fr $(LOCI_LPP_FILES) $(DEPEND_FILES)

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean clean_install distclean spotless),)
-include $(DEPEND_FILES)
endif
