CFLAGS = `gsl-config --cflags`
CFLAGS += -Wall -std=gnu99
LDLIBS = `gsl-config --libs`

all: Vq0.c data plot.tex

Vq0: Vq0.c

data: Vq0
	./$< > $@

plot.tex: data makefile
	echo \
'set terminal epslatex color solid; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$\xi/lt$$"; 		   '\
'set ylabel "$$V(q=0,0)/V_0$$"; 	   '\
'plot [0:][:0] 		   '\
'"$<" using 1:2 notitle with lines lw 4 lc rgb "blue" '\
',"$<" using 1:3 notitle with lines lw 2 lc rgb "red" '\
| tee log | gnuplot



clean:
	$(RM) Vq data plot.tex plot-eps-converted-to.pdf *.eps log
