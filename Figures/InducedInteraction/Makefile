CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: InducedInteraction.tex retneg.tex data log InducedInteraction

InducedInteraction.tex: data makefile
	echo \
'set terminal epslatex color solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$n_F^{1/3}/n_B$$"; 		   '\
'set ylabel "$$2m_F/k_F W_{FF}^\text{ind}(k_F,k_F)$$"; 	   '\
'plot [0:0.8][-2.5:0] 		   '\
'0 notitle with lines lc rgb "black" '\
', "$<" index 1 using ($0.6447):2 notitle with lines dt 2 lw 2 lc rgb "black" '\
',"$<" index 1 using 1:2 notitle with lines lc rgb "black" '\
| gnuplot


retneg.tex: data makefile
	echo \
'set terminal epslatex color solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$n_F^{1/3}/n_B$$"; 		   '\
'set ylabel "$$retneg$$"; 	   '\
'plot [0:0.8][0:100.0] 		   '\
'0 notitle with lines lc rgb "black" '\
', 10.0 notitle with lines lc rgb "red" '\
',"$<" index 1 using 1:3 notitle with lines lc rgb "black" '\
| gnuplot


data log: InducedInteraction
	./$< 1> data 2> log  

InducedInteraction: InducedInteraction.o

clean:
	$(RM) log InducedInteraction InducedInteraction.o *.pdf *.eps *.tex data *.gpi