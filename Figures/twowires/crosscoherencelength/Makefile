CFLAGS = `gsl-config --cflags`
CFLAGS += -Wall -std=gnu99
LDLIBS = `gsl-config --libs`

all: data plot.tex

Vequal: Vequal.c

data: Vequal
	./$< > $@

plot.tex: data makefile
	echo \
'set terminal epslatex color ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k_F\xi$$"; 		   '\
'set ylabel "$$k_Fd_c$$"; 	   '\
'plot  [0:][0:2] 		   '\
' "$<"  using 1:2 notitle with lines lc rgb "black" '\
| tee log | gnuplot


clean:
	$(RM) Vequal data plot.tex plot-eps-converted-to.pdf *.eps log *.gpi