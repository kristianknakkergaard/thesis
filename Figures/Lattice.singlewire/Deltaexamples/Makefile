CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend.tex sitedepend.tex kandsitedepend parameters Delta

kdepend.tex: kandsitedepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'if (!exists("MP_LEFT"))   MP_LEFT = .05; 	'\
'if (!exists("MP_RIGHT"))  MP_RIGHT = .95; 	'\
'if (!exists("MP_BOTTOM")) MP_BOTTOM = .15; 	'\
'if (!exists("MP_TOP"))    MP_TOP = .825; 	'\
'if (!exists("MP_GAP"))    MP_GAP = .05; 	'\
'set  multiplot layout 1,2 rowsfirst margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP; '\
'set grid;       			   		'\
'set tics out;       			    '\
'set xtics -pi, pi/2, pi; 			'\
'set xlabel "$$kd$$";				'\
'set key outside top center Left;			'\
'set arrow from 0,-2.5 to 0,2.5 nohead lc rgb "black"; '\
'set arrow from -pi,0 to pi,0 nohead lc rgb "black"; '\
'set title "$$(i)$$" offset -9, 0;							'\
'set key width -10;									'\
'set ylabel "$$2\\Delta_k, \\varepsilon_k$$";		'\
'plot [-pi:pi][-2.5:2.5]			'\
' "$<" index 0  using 1:(2*$$2) title "$$\Delta_k/t_1 \sim \sin(1kd)$$" with lines lw 2 lc rgb "red" '\
',"$<" index 0  using 1:4 title "$$\\varepsilon_k / t_1$$" with lines lw 2 lc rgb "black";	'\
'set title "$$(ii)$$" offset -9, 0;							'\
'set ytics format " "; 							'\
'unset ylabel; 									'\
'plot [-pi:pi][-2.5:2.5]			'\
' "$<" index 1  using 1:(2*$$2) title "$$\Delta_k/t_1 \sim \sin(1kd)$$" with lines lw 2 dt 2 lc rgb "blue"	'\
',"$<" index 1  using 1:(2*$$3) title "$$\Delta_k/t_1 \sim \sin(2kd)$$" with lines lw 2 lc rgb "blue"		'\
',"$<" index 1  using 1:4 title "$$\\varepsilon_k / t_1$$" with lines lw 2 lc rgb "black";		'\
| gnuplot

sitedepend.tex: kandsitedepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$Delta_k/epsilon_{F,0}$$";'\
'set key bottom right;					'\
'set arrow from 0,-2 to 0,2 nohead lc rgb "black"; '\
'set arrow from 0,0 to 10,0 nohead lc rgb "black"; '\
'plot [0:10][-2:2]			'\
' "$<" index 0  using 4:5 title "$$T/T_F = 0.000$$" with lines lc rgb "#e00000"'\
| gnuplot

kandsitedepend parameters: Delta
	./$< 1> kandsitedepend 2> parameters


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex parameters kandsitedepend *.gpi