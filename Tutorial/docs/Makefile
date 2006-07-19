TEX_FILES = appendix.tex \
		datatype.tex \
		intro.tex \
		io.tex \
		main.tex 


SOURCE_FILES =  entities_cc.tex \
		helpers_cc.tex \
		Makefile_ex.tex

tutorial.pdf: main.dvi
	dvips -Ppdf -o tutorial.ps main
	ps2pdf tutorial.ps
	rm tutorial.ps

main.dvi: ${TEX_FILES} ${SOURCE_FILES}
	latex main; latex main

Makefile_ex.tex: ../1D-Diffusion/Makefile
	cat begin.tex > Makefile_ex.tex
	cat ../1D-Diffusion/Makefile >> Makefile_ex.tex
	echo '\end{verbatim}' >> Makefile_ex.tex

entities_cc.tex: ../Datatypes/entities.cc
	cat begin.tex > entities_cc.tex
	cat ../Datatypes/entities.cc >> entities_cc.tex
	echo '\end{verbatim}' >> entities_cc.tex

helpers_cc.tex: ../Datatypes/helpers.cc
	cat begin.tex > helpers_cc.tex
	cat ../Datatypes/helpers.cc >> helpers_cc.tex
	echo '\end{verbatim}' >> helpers_cc.tex



