
main.dvi: main.tex
	latex main; bibtex main; latex main

Makefile_ex.tex: ../Datatypes/Makefile
	echo "\begin{verbatim}" > Makefile_ex.tex
	cat ../Datatypes/Makefile >> Makefile_ex.tex
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



main.tex : intro.tex dbase.tex containers_cc.tex relations_cc.tex entities_cc.tex helpers_cc.tex Makefile_ex.tex

