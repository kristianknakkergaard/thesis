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
'set xlabel "$$q\xi$$"; 		   '\
'set ylabel "$$V_\text{ind}(q,0)/V_0$$"; 	   '\
'plot [0:0.8][:0] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 1 using 1:2 notitle with lines lc rgb "black" '\
| gnuplot


retneg.tex: data makefile
	echo \
'set terminal epslatex color solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$q\xi$$"; 		   '\
'set ylabel "$$V_\text{ind}(q,0)/V_0$$"; 	   '\
'plot [0:0.8][0:100.0] 		   '\
'0 notitle with lines lc rgb "black" '\
', 20.0 notitle with lines lc rgb "red" '\
',"$<" index 1 using 1:3 notitle with lines lc rgb "black" '\
| gnuplot


data log: InducedInteraction
	./$< 1> data 2> log  

InducedInteraction: InducedInteraction.o

clean:
	$(RM) log InducedInteraction InducedInteraction.o *.pdf *.eps *.tex data *.gpi