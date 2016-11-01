CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: nBdepend.tex nBdepend

nBdepend.tex: nBdepend Makefile
	echo \
'set terminal epslatex solid; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$n_F/n_B^{1/3}$$"; 		   '\
'set ylabel "$$\max_{k}[Delta_k(T=0)]/epsilon_{F,0}$$";   	   '\
'plot [0:0.2][0:1.2] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 1  using 1:2 notitle with lines lc rgb "black"'\
| gnuplot

nBdepend: Delta
	./$< > $@

Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.gpi *.pdf *.eps *depend.tex nBdepend