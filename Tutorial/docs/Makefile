
main.dvi: main.tex
	latex main; bibtex main; latex main

Makefile_ex.tex: ../lfvm/Makefile
	echo "\begin{verbatim}" > Makefile_ex.tex
	cat ../lfvm/Makefile >> Makefile_ex.tex
	echo "\end{verbatim}" >> Makefile_ex.tex


entities_cc.tex: ../Datatypes/entities.cc
	echo "\begin{verbatim}" > entities_cc.tex
	cat ../Datatypes/entities.cc >> entities_cc.tex
	echo "\end{verbatim}" >> entities_cc.tex

containers_cc.tex: ../Datatypes/containers.cc
	echo "\begin{verbatim}" > containers_cc.tex
	cat ../Datatypes/containers.cc >> containers_cc.tex
	echo "\end{verbatim}" >> containers_cc.tex

relations_cc.tex: ../Datatypes/relations.cc
	echo "\begin{verbatim}" > relations_cc.tex
	cat ../Datatypes/relations.cc >> relations_cc.tex
	echo "\end{verbatim}" >> relations_cc.tex

helpers_cc.tex: ../Datatypes/helpers.cc
	echo "\begin{verbatim}" > helpers_cc.tex
	cat ../Datatypes/helpers.cc >> helpers_cc.tex
	echo "\end{verbatim}" >> helpers_cc.tex

grid_reader_cc.tex: ../lfvm/grid_reader.cc
	echo "\begin{verbatim}" > grid_reader_cc.tex
	cat ../lfvm/grid_reader.cc >> grid_reader_cc.tex
	echo "\end{verbatim}" >> grid_reader_cc.tex

setup_edges_cc.tex: ../lfvm/setup_edges.cc
	echo "\begin{verbatim}" > setup_edges_cc.tex
	cat ../lfvm/setup_edges.cc >> setup_edges_cc.tex
	echo "\end{verbatim}" >> setup_edges_cc.tex

cell_props_cc.tex: ../lfvm/cell_props.cc
	echo "\begin{verbatim}" > cell_props_cc.tex
	cat ../lfvm/cell_props.cc >> cell_props_cc.tex
	echo "\end{verbatim}" >> cell_props_cc.tex

area_reduce_cc.tex: ../lfvm/area_reduce.cc
	echo "\begin{verbatim}" > area_reduce_cc.tex
	cat ../lfvm/area_reduce.cc >> area_reduce_cc.tex
	echo "\end{verbatim}" >> area_reduce_cc.tex



main.dvi : intro.tex dbase.tex containers_cc.tex relations_cc.tex entities_cc.tex helpers_cc.tex Makefile_ex.tex grid_reader_cc.tex setup_edges_cc.tex cell_props_cc.tex area_reduce_cc.tex




