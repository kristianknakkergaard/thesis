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
'set ylabel "$$(\mu + t_2)/t_1$$";				'\
'set key top left;					'\
'set arrow from 0,0 to 1,0 nohead lc rgb "black"; '\
'plot [0:1][-2:3]			'\
' "$<" index 1  using 1:2 title "$$t_2 / t_1 = 0.00$$" with lines lw 2 lc rgb "black"'\
',"$<" index 1  using (1 -$$1):(-$$2) notitle with lines lw 2 lc rgb "black"'\
',"$<" index 2  using 1:(0.25 + $$2) title "$$t_2 / t_1 = 0.25$$" with lines lw 2 lc rgb "#e00000"'\
',"$<" index 3  using 1:(0.5 + $$2) title "$$t_2 / t_1 = 0.50$$" with lines lw 2 lc rgb "#167000"'\
',"$<" index 4  using 1:(1.0 + $$2) title "$$t_2 / t_1 = 1.00$$" with lines lw 2 lc rgb "#9100aa"'\
| gnuplot

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex ndepend *.gpi

	