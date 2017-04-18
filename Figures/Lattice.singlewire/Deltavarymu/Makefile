CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: mudepend.tex mudepend Edepend.tex parameters 

mudepend.tex: mudepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$t_2 / t_1$$"; 		    '\
'set ylabel "$$CS_1, E_{\min}$$";'\
'set key top left;					'\
'set arrow from 0,-2.1 to 0,2.1 nohead lc rgb "black"; '\
'set arrow from -1,0 to 2,0 nohead lc rgb "black"; '\
'plot [-1:2][-2.1:2.1]			'\
' "$<" index 0  using 1:2 notitle with lines lc rgb "red"'\
',"$<" index 0  using 1:5 notitle with lines lc rgb "blue"'\
| gnuplot

Edepend.tex: mudepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$t_2 / t_1$$"; 		    '\
'set ylabel "$$E_{\min}$$";'\
'set key top left;					'\
'set arrow from 0,-2.1 to 0,2.1 nohead lc rgb "black"; '\
'set arrow from -1,0 to 3,0 nohead lc rgb "black"; '\
'plot [-1:2][-2.1:2.1]			'\
' "$<" index 0  using 1:5 notitle with lines lc rgb "red"'\
| gnuplot

mudepend parameters: Delta
	./$< 1> mudepend 2> parameters


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex parameters kdepend *.gpi