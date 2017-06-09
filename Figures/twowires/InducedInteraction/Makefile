CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: InducedInteraction.tex data log InducedInteraction

InducedInteraction.tex: data makefile
	echo \
'set terminal epslatex color solid ; 	'\
'set output "$@"; 			   			'\
'set grid;       			   			'\
'set tics out;       			   		'\
'set key top right Left; 				'\
'set key width -2; 						'\
'set key opaque; 						'\
'set xlabel "$$n_F^{1/3}/n_B$$"; 		'\
'set ylabel "$$2m_F/k_F W_{FF}^\text{ind}(k_F,k_F)$$"; 	   '\
'set arrow from 0.6447,-2.5 to 0.6447,0 nohead dt 2 lw 2 lc rgb "black"; '\
'plot [0:1.5][-2.5:0] 		   '\
'0 notitle with lines lc rgb "black" '\
', "$<" index 0 using 1:2 title "$$Intrawire$$" with lines lc rgb "red" '\
', "$<" index 0 using 1:3 title "$$Interwire$$" with lines lc rgb "blue" '\
| gnuplot

data log: InducedInteraction
	./$< 1> data 2> log  

InducedInteraction: InducedInteraction.o

clean:
	$(RM) log InducedInteraction InducedInteraction.o *.pdf *.eps *.tex data *.gpi