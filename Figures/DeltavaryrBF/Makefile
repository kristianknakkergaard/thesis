CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: rBFdepend.tex rBFdepend 

rBFdepend.tex: rBFdepend Makefile
	echo \
'set terminal epslatex solid; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k_F a_B$$"; 		   '\
'set ylabel "$$\max_{k}[Delta_k(T=0)]/epsilon_{F,0}$$";   	   '\
'plot [0.03:0.04][0:7] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 1  using 1:2 notitle with lines lc rgb "black"'\
| gnuplot

rBFdepend: Delta
	./$< > $@

Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.gpi *.pdf *.eps *depend.tex rBFdepend