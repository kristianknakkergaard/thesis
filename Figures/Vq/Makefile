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
'set xlabel "$$q/k_F$$"; 		   '\
'set ylabel "$$2m_F/k_F V^\text{ind}_{FF}(q,0)$$"; 	   '\
'set key bottom left;					'\
'set arrow from 0,0 to 0,-20 nohead lc rgb "black"; '\
'plot [-20:20][-20:0] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 3 using 2:3 title "$$k_Fl_t = 0.20$$" with lines lc rgb "#002daa" '\
',"$<" index 2 using 2:3 title "$$k_Fl_t = 0.15$$" with lines lc rgb "#167000" '\
',"$<" index 1 using 2:3 title "$$k_Fl_t = 0.10$$" with lines lc rgb "#e05900" '\
',"$<" index 0 using 2:3 title "$$k_Fl_t = 0.05$$" with lines lc rgb "#e00000" '\
| tee log | gnuplot



clean:
	$(RM) Vq data plot.tex plot-eps-converted-to.pdf *.eps log
