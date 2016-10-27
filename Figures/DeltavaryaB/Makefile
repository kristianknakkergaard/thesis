CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: aBdepend.tex mudepend.tex aBdepend 

aBdepend.tex: aBdepend Makefile
	echo \
'set terminal epslatex solid; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k_F a_B$$"; 		   '\
'set ylabel "$$\max_{k}[Delta_k(T=0)]/epsilon_{F,0}$$";   	   '\
'plot [0:0.3][0:1.2] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 0  using 1:2 notitle with lines lc rgb "black"'\
| gnuplot

mudepend.tex: aBdepend Makefile
	echo \
'set terminal epslatex solid; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k_F a_B$$"; 		   '\
'set ylabel "$$\epsilon_F/epsilon_{F,0}$$";   	   '\
'plot [0:0.3][0.936:0.944] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 0  using 1:3 notitle with lines lc rgb "black"'\
| gnuplot



aBdepend: Delta
	./$< > $@


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.gpi *.pdf *.eps *depend.tex aBdepend