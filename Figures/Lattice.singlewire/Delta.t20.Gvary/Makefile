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
'plot [-pi:pi][-3:3]			'\
' "$<" index 0  using 1:2 title "$$G = 0.5$$" with lines lc rgb "#e00000"'\
',"$<" index 1  using 1:2 title "$$G = 1.0$$" with lines lc rgb "#e05900"'\
',"$<" index 2  using 1:2 title "$$G = 2.0$$" with lines lc rgb "#167000"'\
',"$<" index 3  using 1:2 title "$$G = 4.0$$" with lines lc rgb "#002daa"'\
',"$<" index 4  using 1:2 title "$$G = 8.0$$" with lines lc rgb "#9100aa"'\
',"$<" index 4  using 1:(0.92 * sin($$1)) notitle with lines lw 2 lc rgb "black"'\
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
' "$<" index 0  using 4:5 title "$$G = 0.5$$" with lines lc rgb "#e00000"'\
',"$<" index 1  using 4:5 title "$$G = 1.0$$" with lines lc rgb "#e05900"'\
',"$<" index 2  using 4:5 title "$$G = 2.0$$" with lines lc rgb "#167000"'\
',"$<" index 3  using 4:5 title "$$G = 4.0$$" with lines lc rgb "#002daa"'\
',"$<" index 4  using 4:5 title "$$G = 8.0$$" with lines lc rgb "#9100aa"'\
| gnuplot

kandsitedepend parameters: Delta
	./$< 1> kandsitedepend 2> parameters


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex parameters kandsitedepend *.gpi