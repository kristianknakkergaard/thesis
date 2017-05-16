CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend.tex kdepend parameters 

kdepend.tex: kdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$Delta_k/epsilon_{F,0}$$";'\
'set key top left;					'\
'set key width -1;					'\
'set arrow from 0,-0.3 to 0,0.3 nohead lc rgb "black"; '\
'set arrow from -10,0 to 10,0 nohead lc rgb "black"; '\
'plot [-10:10][-0.3:0.3]			'\
' "$<" index 0  using 1:2 title "$$T/T_F = 0.000$$" with lines lc rgb "#e00000"'\
',"$<" index 1  using 1:2 title "$$T/T_F = 0.075$$" with lines lc rgb "#e05900"'\
',"$<" index 2  using 1:2 title "$$T/T_F = 0.100$$" with lines lc rgb "#167000"'\
',"$<" index 3  using 1:2 title "$$T/T_F = 0.120$$" with lines lc rgb "#002daa"'\
',"$<" index 4  using 1:2 title "$$T/T_F = 0.127$$" with lines lc rgb "#9100aa"'\
| gnuplot

kdepend parameters: Delta
	./$< 1> kdepend 2> parameters


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex parameters kdepend *.gpi