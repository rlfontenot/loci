#LOCI_BASE = /usr/tmp/qxue/Loci
#LOCI_BASE ?= /simcenter/data4/loci/opt/Loci-Linux-i686-mpic++-rel-3-1-exp-7
#LOCI_BASE ?= /simcenter/data4/loci/opt/Loci-Linux-i686-mpic++-08.09.07
include $(LOCI_BASE)/Loci.conf

#COPTLESS =	-O0 -g -ansi -pedantic
COPTLESS =	-O0 -g 
COPT =		$(COPTLESS)

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
		get_hex_fine_grid.o \
		get_prism_fine_grid.o \
		fine_grid_output.o  \
		merge_mixed_face.o \
		transfer_plan.o \
		get_mixed_face_fine_grid.o \
		refmesh.o 




MODULES =	

LOCAL_LIBS =	

JUNK =		*~ \
		core \
		ti_files \
		ii_files \
		rii_files

INCLUDES=	-I./include

default: mymarker myrefmesh 

mymarker: $(MARKOBJS) $(MODULES) 
	$(LD) -o  mymarker  $(MARKOBJS) $(MODULES) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

myrefmesh:$(REFOBJS) $(MODULES) 
	$(LD) -o  myrefmesh  $(REFOBJS) $(MODULES) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)


FRC :

clean:
	rm -fr $(MARKOBJS) $(REFOBJS)  $(TARGET) $(JUNK)

distclean: clean
	rm -fr $(DEPEND_FILES) myrefmesh mymarker

DEPEND_FILES=$(subst .o,.d,$(MARKOBJS) $(REFOBJS))
#include automatically generated dependencies
include $(DEPEND_FILES)
