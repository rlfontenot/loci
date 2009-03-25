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
include $(LOCI_BASE)/Loci.conf

#COPTLESS =	-O0 -g -ansi -pedantic
#COPT =		$(COPTLESS)

MARKOBJS =	classify_cell.o \
		transfer_fc.o \
		quadface.o  \
		face.o \
		build_prismcell.o \
		build_hexcell.o \
		build_general_cell.o \
		create_prism_maps.o \
		create_hexcell_maps.o \
		extract_hex_edge.o \
		extract_hex_face.o \
		extract_prism_face.o \
		diamondcell.o \
		prism.o \
		hexcell.o \
		merge_general_face.o \
		merge_general_edge.o \
		merge_hex_face.o \
		merge_prism_face.o \
		set_hex_nums.o \
		set_general_nums.o \
		set_prism_nums.o \
		globals.o \
		node_edge.o \
		make_general_cellplan.o \
		make_hex_cellplan.o \
		make_prism_cellplan.o \
		update_quad_face.o \
		update_general_face.o \
		update_prism_face.o \
		update_general_edge.o \
		balance_hex_cell.o \
		balance_prism_cell.o \
		balance_general_cell.o \
		cellplan_io.o \
		check_folded_face.o \
		get_node_remap.o \
		read_tag.o \
		merge_mixed_face.o \
		update_mixed_face.o \
		transfer_plan.o \
		parse_tree.o \
		marker.o

REFOBJS =	classify_cell.o \
		transfer_fc.o \
		quadface.o \
		face.o \
		build_hexcell.o \
		build_prismcell.o \
		build_general_cell.o \
		create_hexcell_maps.o \
		create_prism_maps.o \
		extract_hex_edge.o \
		extract_hex_face.o \
		extract_prism_face.o \
		diamondcell.o \
		hexcell.o \
		prism.o \
		merge_general_face.o \
		merge_general_edge.o \
		merge_hex_face.o \
		merge_prism_face.o \
		set_hex_nums.o \
		set_general_nums.o \
		set_prism_nums.o \
		get_c1_general.o \
		get_c1_hex.o \
		get_c1_prism.o \
		globals.o \
		node_edge.o \
		get_node_remap.o \
		cellplan_io.o \
		set_offset.o \
		get_fine_grid.o \
		get_fine_face.o \
		get_hex_fine_grid.o \
		get_hex_fine_face.o \
		get_prism_fine_grid.o \
		get_prism_fine_face.o \
		write_vog.o  \
		color_matrix.o \
		merge_mixed_face.o \
		transfer_plan.o \
		get_mixed_face_fine_grid.o \
		refmesh.o 




#MODULES =	-lxml2	

#LOCAL_LIBS =	

LIBXML2_LIBS ?= -lxml2
LIBXML2_DEFINE ?= -DUSE_LIBXML2
LIBXML2_INCLUDE ?= -I/usr/include/libxml2

LOCAL_LIBS      = $(LIBXML2_LIBS)

JUNK =		*~ \
		core \
		ti_files \
		ii_files \
		rii_files

INCLUDES=	-I./include $(LIBXML2_INCLUDE) $(LIBXML2_DEFINE)

default: marker refmesh 

marker: $(MARKOBJS) $(MODULES) 
	$(LD) -o  marker  $(MARKOBJS) $(MODULES) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

refmesh:$(REFOBJS) $(MODULES) marker
	$(LD) -o  refmesh  $(REFOBJS) $(MODULES) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)


FRC :

clean:
	rm -fr $(MARKOBJS) $(REFOBJS)  $(TARGET) $(JUNK)

distclean: clean
	rm -fr $(DEPEND_FILES) refmesh marker

DEPEND_FILES=$(subst .o,.d,$(MARKOBJS) $(REFOBJS))
#include automatically generated dependencies
include $(DEPEND_FILES)
