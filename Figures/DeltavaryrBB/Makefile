CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: rBBdepend.tex rBBdepend 

rBBdepend.tex: rBBdepend Makefile
	echo \
'set terminal epslatex solid; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k_F a_B$$"; 		   '\
'set ylabel "$$\max_{k}[Delta_k(T=0)]/epsilon_{F,0}$$";   	   '\
'plot [0:0.01][0:30] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 1  using 1:2 notitle with lines lc rgb "black"'\
| gnuplot

rBBdepend: Delta
	./$< > $@

Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.gpi *.pdf *.eps *depend.tex rBBdepend