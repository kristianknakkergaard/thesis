CFLAGS = `gsl-config --cflags`
CFLAGS += -Wall -std=gnu99
LDLIBS = `gsl-config --libs`

all: Vx.c data plot.tex

Vx: Vx.c

data: Vx
	./$< > $@

plot.tex: data makefile
	echo \
'set terminal epslatex color solid ; 								'\
'set output "$@"; 			   										'\
'set grid;       			   										'\
'set tics out;       			  				 					'\
'set xlabel "$$k_Fx$$"; 		   									'\
'set ylabel "$$\tilde{V}_{FF}^\text(x,0)/\epsilon_{F,0}$$"; 	   	'\
'set key bottom left;												'\
'set key width -2;													'\
'set arrow from 0,0 to 0,-60 nohead lc rgb "black"; 				'\
'plot [-5:5][-3:0] 		   										'\
' "$<" index 0 using 2:4 notitle with lines lw 2 lc rgb "black" 				'\
',"$<" index 2 using 2:3 title "$$k_Fl_t = 0.6$$" with lines lc rgb "#002daa" 	'\
',"$<" index 1 using 2:3 title "$$k_Fl_t = 0.4$$" with lines lc rgb "#167000" 	'\
',"$<" index 0 using 2:3 title "$$k_Fl_t = 0.2$$" with lines lc rgb "#e00000" 	'\
| tee log | gnuplot



clean:
	$(RM) Vx data plot.tex plot-eps-converted-to.pdf *.eps log
