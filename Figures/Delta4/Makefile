CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: Tdepend.tex Tdepend

Tdepend.tex: Tdepend Makefile
	echo \
'set terminal epslatex solid; 							'\
'set output "$@"; 			   							'\
'set grid;       			   							'\
'set tics out;       			   						'\
'set key top right; 									'\
'set key width -3;										'\
'set xlabel "$$T/T_F$$"; 		   						'\
'set ylabel "$$\max_k[\Delta_k]/\epsilon_{F,0}$$"; 	   	'\
'plot [0:0.18][0:0.4] 		   							'\
' "$<" index 1 using 1:3 title "$$Numerical data$$" with lines lc rgb "black" 													'\
',"$<" index 1 using 1:(1.6789 * 0.279949 * sqrt(1.0 - $$1/0.1294)) title "$$Asymptotic behaviour$$" with lines lc rgb "blue" 	'\
',"$<" index 2 using 1:2 title "$$Approximate behaviour$$" with lines lc rgb "red" 												'\
| gnuplot


clean:
	$(RM) Tdepend.tex