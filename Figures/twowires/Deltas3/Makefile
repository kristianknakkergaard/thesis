CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: nBdepend.tex nBdepend

nBdepend.tex: nBdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$n_F/n_B^{1/3}$$"; 		    '\
'set ylabel "$$d_equal$$";'\
'set key top left;					'\
'plot [0.01:0.65][0:1]			'\
' "$<" index 0  using 1:2 notitle with lines lc rgb "black"’\
| gnuplot

nBdepend: Delta
	./$< 1> nBdepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex Tdepend kdepend *.gpi