CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend.tex kdepend parameters Delta

kdepend.tex: kdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$Delta_k/epsilon_{F,0}$$";'\
'set key top left;					'\
'set arrow from 0,-0.5 to 0,0.5 nohead lc rgb "black"; '\
'set arrow from -3.2,0 to 3.2,0 nohead lc rgb "black"; '\
'plot [-3.14:3.14][-0.5:0.5]			'\
' "$<" index 0  using 1:2 title "$$T/T_F = 0.000$$" with lines lc rgb "#e00000"'\
',"$<" index 0  using 1:(0.305 * sin($$1) + 0.0 * sin(2*$$1)) notitle with lines lc rgb "black"'\
| gnuplot

kdepend parameters: Delta
	./$< 1> kdepend 2> parameters


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex parameters kdepend *.gpi