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
'set xlabel "$$1/k_Fl_t$$"; 		   '\
'set ylabel "$$V_{FF}^\text{ind}(q=0,0)/V_0$$"; 	   '\
'plot [0:][:0] 		   '\
'"$<" using 1:2 notitle with lines lw 4 lc rgb "blue" '\
',"$<" using 1:3 notitle with lines lw 2 lc rgb "red" '\
| tee log | gnuplot



clean:
	$(RM) Vq0 data plot.tex plot-eps-converted-to.pdf *.eps *.o *.log
