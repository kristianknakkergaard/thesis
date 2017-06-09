CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend dispersion kdepend.tex

kdepend.tex: kdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'if (!exists("MP_LEFT"))   MP_LEFT = .05; 	'\
'if (!exists("MP_RIGHT"))  MP_RIGHT = .95; 	'\
'if (!exists("MP_BOTTOM")) MP_BOTTOM = .15; 	'\
'if (!exists("MP_TOP"))    MP_TOP = .9; 	'\
'if (!exists("MP_GAP"))    MP_GAP = .05; 	'\
'set  multiplot layout 2,2 rowsfirst margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP; '\
'unset grid;       			   												'\
'set tics out;       			   											'\
'set xtics -pi, pi/2, pi; 													'\
'set ytics -3, 1.5, 3; 														'\
'set arrow from 0,-3 to 0,3 nohead lc rgb "black"; 							'\
'set arrow from -pi,0 to pi,0 nohead lc rgb "black"; 				'\
'set ylabel "$$\varepsilon_k, \Delta_k$$";									'\
'set xtics format " ";														'\
'unset xlabel;																'\
'set title "$$nu = 1$$" offset -8.0, -2.5;							'\
'plot  [-pi:pi][-3:3]												'\
' "$<"  using 1:2 notitle with lines lc rgb "blue"							'\
',"$<"  using 1:4 notitle with lines lc rgb "red";							'\
'unset ylabel;																'\
'set ytics format " ";														'\
'set title "$$nu = 1$$"; 											'\
'plot  [-pi:pi][-3:3]												'\
' "$<"  using 1:2 notitle with lines lc rgb "blue"							'\
',"$<"  using 1:5 notitle with lines lc rgb "red";							'\
'set ylabel "$$\varepsilon_k, \Delta_k$$";									'\
'set xlabel "$$kd$$"; 		   								 				'\
'set title "$$nu = 0$$"; 											'\
'plot  [-pi:pi][-3:3]												'\
' "$<"  using 1:3 notitle with lines lc rgb "blue"							'\
',"$<"  using 1:4 notitle with lines lc rgb "red";							'\
'unset ylabel;																'\
'set ytics format " ";														'\
'set title "$$nu = 2$$"; 											'\
'plot  [-pi:pi][-3:3]												'\
' "$<"  using 1:3 notitle with lines lc rgb "blue"							'\
',"$<"  using 1:5 notitle with lines lc rgb "red";     						'\
'unset multiplot; 															'\
| gnuplot

kdepend: dispersion
	./$< 1> kdepend

dispersion: dispersion.o

clean:
	$(RM) dispersion dispersion.o *.pdf *.eps *depend.tex kdepend *.gpi


 