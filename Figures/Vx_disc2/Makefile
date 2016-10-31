CFLAGS = `gsl-config --cflags`
CFLAGS += -Wall -std=gnu99
LDLIBS = `gsl-config --libs`

all: data plot.tex

Disc: Disc.c

data: Disc
	./$< > $@

plot.tex: data makefile
	echo \
'set terminal epslatex color solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k_Fl_t$$"; 		   '\
'set ylabel "$$k_Fd$$"; 	   '\
'plot [0:][0:] 		   '\
' "$<" using 1:2 notitle with lines lw 2 lc rgb "blue" '\
',"$<" using 1:4 notitle with lines lw 2 lc rgb "red" '\
| tee log | gnuplot


clean:
	$(RM) Disc data plot.tex plot-eps-converted-to.pdf *.eps log *.gpi