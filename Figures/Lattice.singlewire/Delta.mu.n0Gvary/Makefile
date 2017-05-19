CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: ndepend.tex

ndepend.tex: ndepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$n$$"; 		  		'\
'set ylabel "$$\mu$$";				'\
'set key top left;					'\
'set arrow from 0,0 to 1,0 nohead lc rgb "black"; '\
'plot [0:1][-2.5:2.5]			'\
' "$<" index 0  using 1:(-cos(3.14159265 * $$1)) title "$$G = 0.0$$" with lines lw 2 lc rgb "black"'\
',"$<" index 2  using 1:2 title "$$G = 2.0$$" with lines lw 2 lc rgb "#e00000"'\
',"$<" index 2  using (1 -$$1):(-$$2) notitle with lines lw 2 lc rgb "#e00000"'\
',"$<" index 4  using 1:2 title "$$G = 4.0$$" with lines lw 2 lc rgb "#167000"'\
',"$<" index 4  using (1 -$$1):(-$$2) notitle with lines lw 2 lc rgb "#167000"'\
',"$<" index 6  using 1:2 title "$$G = 8.0$$" with lines lw 2 lc rgb "#9100aa"'\
',"$<" index 6  using (1 -$$1):(-$$2) notitle with lines lw 2 lc rgb "#9100aa"'\
| gnuplot

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex ndepend *.gpi