TEX_FILES	=	appendix.tex \
				datatype.tex \
				intro.tex    \
				io.tex       \
				main.tex


SOURCE_FILES	=	entities_cpp.tex \
					helpers_cpp.tex \
					Makefile_ex.tex

tutorial.pdf: ${TEX_FILES} ${SOURCE_FILES}
	pdflatex  main && pdflatex main && mv main.pdf tutorial.pdf


clean:
	./clean

Makefile_ex.tex: ../heat/Makefile
	cat begin.tex        >  Makefile_ex.tex
	cat ../heat/Makefile >> Makefile_ex.tex
	cat end.tex          >> Makefile_ex.tex

entities_cpp.tex: ../Datatypes/entities.cpp
	cat begin.tex                 >  entities_cpp.tex
	cat ../Datatypes/entities.cpp >> entities_cpp.tex
	cat end.tex                   >> entities_cpp.tex

helpers_cpp.tex: ../Datatypes/helpers.cpp
	cat begin.tex                >  helpers_cpp.tex
	cat ../Datatypes/helpers.cpp >> helpers_cpp.tex
	cat end.tex                  >> helpers_cpp.tex



