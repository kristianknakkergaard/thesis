CFLAGS = `gsl-config --cflags`
CFLAGS += -Wall -std=gnu99
LDLIBS = `gsl-config --libs`

all: Vq.c data plot.tex

Vq: Vq.c

data: Vq
	./$< > $@

plot.tex: data makefile
	echo \
'set terminal epslatex color solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$q\xi$$"; 		   '\
'set ylabel "$$V_\text{ind}(q,0)/V_0$$"; 	   '\
'plot [-8:8][-4.3:0] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 0 using 2:3 notitle with lines lw 2 lc rgb "red" '\
',"$<" index 1 using 2:3 notitle with lines lw 2 lc rgb "blue"  '\
',"$<" index 2 using 2:3 notitle with lines lw 2 lc rgb "red" '\
',"$<" index 3 using 2:3 notitle with lines lw 2 lc rgb "blue"  '\
',"$<" index 4 using 2:3 notitle with lines lw 2 lc rgb "red" '\
',"$<" index 5 using 2:3 notitle with lines lw 2 lc rgb "blue"  '\
',"$<" index 6 using 2:3 notitle with lines lw 2 lc rgb "red" '\
| tee log | gnuplot



clean:
	$(RM) Vq data plot.tex plot-eps-converted-to.pdf *.ps log
