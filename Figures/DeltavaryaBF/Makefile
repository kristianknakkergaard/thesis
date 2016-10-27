CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: aBFdepend.tex aBFdepend 

aBFdepend.tex: aBFdepend Makefile
	echo \
'set terminal epslatex solid; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k_F a_{BF}$$"; 		   '\
'set ylabel "$$\max_{k}[\Delta_k(T=0)]/\epsilon_{F,0}$$";   	   '\
'set key top left;		   	   '\
'plot [0.15:0.2][0:0.25] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 0  using 1:2 notitle with lines lc rgb "black"'\
| gnuplot

aBFdepend: Delta
	./$< > $@


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.gpi *.pdf *.eps *depend.tex aBFdepend