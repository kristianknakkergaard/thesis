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
'set ylabel "$$\tilde{V}^{11}_{\text{ind}}(x,0)}/\epsilon_{F,0}$$"; 	   	'\
'set key bottom left;												'\
'set key width -2;													'\
'set arrow from 0,0 to 0,-30 nohead lc rgb "black"; 				'\
'plot [-1:1][-30:0] 		   										'\
' "$<" index 6 using 2:3 title "$$k_Fd = 0.3$$" with lines 		lc rgb "#002daa" '\
',"$<" index 4 using 2:3 title "$$k_Fd = 0.2$$" with lines 		lc rgb "#167000" '\
',"$<" index 2 using 2:3 title "$$k_Fd = 0.1$$" with lines 		lc rgb "#e00000" '\
',"$<" index 0 using 2:3 title "$$k_Fd = 0.0$$" with lines lw 2 lc rgb "black" '\
| tee log | gnuplot



clean:
	$(RM) Vx data plot.tex plot-eps-converted-to.pdf *.eps log
