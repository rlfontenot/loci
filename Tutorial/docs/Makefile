TEX_FILES =	appe.tex \
		appendix.tex \
		basic.tex \
		basic1.tex \
		bv.tex \
		datatype.tex \
		dbase.tex \
		example.tex \
		heat.tex \
		intro.tex \
		io.tex \
		lap.tex \
		lists.tex \
		main.tex \
		mesh.tex \
		pi.tex \
		pi1.tex \
		rep.tex \
		start.tex \
		tech.tex \
		title.tex

IMAGE_FILES =	map4.eps \
		pi.eps \
		figures/edge_area.eps \
		figures/nested.eps

SOURCE_FILES =	area_reduce_cc.tex \
		cell_props_cc.tex \
		containers_cc.tex \
		entities_cc.tex \
		gradient_cc.tex \
		grid_reader_cc.tex \
		helpers_cc.tex \
		integration_cc.tex \
		Makefile_ex.tex \
		petsc1_main_cc.tex \
		petsc1_rules_cc.tex \
		relations_cc.tex \
		setup_edges_cc.tex

main.dvi: ${TEX_FILES} ${IMAGE_FILES} ${SOURCE_FILES}
	latex main; bibtex main; latex main

Makefile_ex.tex: ../lfvm/Makefile
	cat bv.tex > Makefile_ex.tex
	cat ../lfvm/Makefile >> Makefile_ex.tex
	echo '\end{verbatim}' >> Makefile_ex.tex

entities_cc.tex: ../Datatypes/entities.cc
	cat bv.tex > entities_cc.tex
	cat ../Datatypes/entities.cc >> entities_cc.tex
	echo '\end{verbatim}' >> entities_cc.tex

containers_cc.tex: ../Datatypes/containers.cc
	cat bv.tex > containers_cc.tex
	cat ../Datatypes/containers.cc >> containers_cc.tex
	echo '\end{verbatim}' >> containers_cc.tex

relations_cc.tex: ../Datatypes/relations.cc
	cat bv.tex > relations_cc.tex
	cat ../Datatypes/relations.cc >> relations_cc.tex
	echo '\end{verbatim}' >> relations_cc.tex

helpers_cc.tex: ../Datatypes/helpers.cc
	cat bv.tex > helpers_cc.tex
	cat ../Datatypes/helpers.cc >> helpers_cc.tex
	echo '\end{verbatim}' >> helpers_cc.tex

grid_reader_cc.tex: ../lfvm/grid_reader.cc
	cat bv.tex > grid_reader_cc.tex
	cat ../lfvm/grid_reader.cc >> grid_reader_cc.tex
	echo '\end{verbatim}' >> grid_reader_cc.tex

setup_edges_cc.tex: ../lfvm/setup_edges.cc
	cat bv.tex > setup_edges_cc.tex
	cat ../lfvm/setup_edges.cc >> setup_edges_cc.tex
	echo '\end{verbatim}' >> setup_edges_cc.tex

cell_props_cc.tex: ../lfvm/cell_props.cc
	cat bv.tex > cell_props_cc.tex
	cat ../lfvm/cell_props.cc >> cell_props_cc.tex
	echo '\end{verbatim}' >> cell_props_cc.tex

area_reduce_cc.tex: ../lfvm/area_reduce.cc
	cat bv.tex > area_reduce_cc.tex
	cat ../lfvm/area_reduce.cc >> area_reduce_cc.tex
	echo '\end{verbatim}' >> area_reduce_cc.tex

integration_cc.tex: ../lfvm/integration.cc
	cat bv.tex > integration_cc.tex
	cat ../lfvm/integration.cc >> integration_cc.tex
	echo '\end{verbatim}' >> integration_cc.tex

gradient_cc.tex: ../lfvm/gradient.cc
	cat bv.tex > gradient_cc.tex
	cat ../lfvm/gradient.cc >> gradient_cc.tex
	echo '\end{verbatim}' >> gradient_cc.tex

petsc1_main_cc.tex: ../petsc/petsc1_main.cc
	cat bv.tex > petsc1_main_cc.tex
	cat ../petsc/petsc1_main.cc >> petsc1_main_cc.tex
	echo '\end{verbatim}' >> petsc1_main_cc.tex

petsc1_rules_cc.tex: ../petsc/petsc1_rules.cc
	cat bv.tex > petsc1_rules_cc.tex
	cat ../petsc/petsc1_rules.cc >> petsc1_rules_cc.tex
	echo '\end{verbatim}' >> petsc1_rules_cc.tex


