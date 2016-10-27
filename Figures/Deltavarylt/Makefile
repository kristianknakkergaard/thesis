CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: ltdepend.tex ltdepend 

ltdepend.tex: ltdepend Makefile
	echo \
'set terminal epslatex; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k_Fl_t$$"; 		   '\
'set ylabel "$$\max_{k}[\Delta_k(T=0)]/\epsilon_{F,0}$$";   	   '\
'set key top left;		   	   '\
'plot [0:0.3][0:0.31] 		   '\
' "$<" index 0  using 1:2 notitle with lines lc rgb "black"'\
| gnuplot

ltdepend: Delta
	./$< > $@


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.gpi *.pdf *.eps *depend.tex ltdepend