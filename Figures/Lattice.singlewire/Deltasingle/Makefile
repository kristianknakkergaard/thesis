CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend.tex sitedepend.tex kandsitedepend parameters Delta

kdepend.tex: kandsitedepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$Delta_k/epsilon_{F,0}$$";'\
'set key bottom right;					'\
'set arrow from 0,-3 to 0,3 nohead lc rgb "black"; '\
'set arrow from -3.2,0 to 3.2,0 nohead lc rgb "black"; '\
'plot [-3.14:3.14][-3:3]			'\
' "$<" index 0  using 1:2 title "$$T/T_F = 0.000$$" with lines lc rgb "#e00000"'\
',"$<" index 0  using 1:3 title "$$\varepsilon_k / t_1$$" with lines lc rgb "black"'\
| gnuplot

sitedepend.tex: kandsitedepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$Delta_k/epsilon_{F,0}$$";'\
'set key bottom right;					'\
'set arrow from 0,-2 to 0,2 nohead lc rgb "black"; '\
'set arrow from 0,0 to 10,0 nohead lc rgb "black"; '\
'plot [0:10][-2:2]			'\
' "$<" index 0  using 4:5 title "$$T/T_F = 0.000$$" with lines lc rgb "#e00000"'\
| gnuplot

kandsitedepend parameters: Delta
	./$< 1> kandsitedepend 2> parameters


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex parameters kandsitedepend *.gpi