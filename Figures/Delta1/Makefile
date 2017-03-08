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
'plot [-10:10][-5:5]			'\
' "$<" index 0     using 1:2 title "$$T/T_F = 0.000$$" with lines lc rgb "#e00000"'\
| gnuplot

kdepend parameters: Delta
	./$< 1> kdepend 2> parameters


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex parameters kdepend *.gpi