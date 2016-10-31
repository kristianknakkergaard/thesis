CFLAGS = `gsl-config --cflags`
CFLAGS += -Wall -std=gnu99
LDLIBS = `gsl-config --libs`

all: Vx.c data plot.tex

Vx: Vx.c

data: Vx
	./$< > $@

plot.tex: data makefile
	echo \
'set terminal epslatex color solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k_Fx$$"; 		   '\
'set ylabel "$$\tilde{V}_{FF}^\text(x,0)/V_0$$"; 	   '\
'plot [-1:1][-30:0] 		   '\
' "$<" index 0 using 2:4 notitle with lines lw 2 lc rgb "black" '\
',"$<" index 0 using 2:3 notitle with lines lw 2 lc rgb "red" '\
',"$<" index 1 using 2:3 notitle with lines lw 2 lc rgb "blue"  '\
',"$<" index 2 using 2:3 notitle with lines lw 2 lc rgb "red" '\
',"$<" index 3 using 2:3 notitle with lines lw 2 lc rgb "blue"  '\
',"$<" index 4 using 2:3 notitle with lines lw 2 lc rgb "red" '\
',"$<" index 5 using 2:3 notitle with lines lw 2 lc rgb "blue"  '\
',"$<" index 6 using 2:3 notitle with lines lw 2 lc rgb "red" '\
| tee log | gnuplot



clean:
	$(RM) Vx data plot.tex plot-eps-converted-to.pdf *.eps log
