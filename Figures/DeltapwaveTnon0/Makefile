CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend.tex Tdepend.tex kdepend Tdepend

kdepend.tex: kdepend Makefile
	echo \
'set terminal epslatex; 	   '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       		   '\
'set key top left;        '\
'set xlabel "$$k/k_F$$"; 		   '\
'set ylabel "$$\Delta_k/\epsilon_F$$"; 	   '\
'plot [-10:10][-0.3:0.3] 		   '\
' "$<" index 1  using 1:2 title "$$T/T_F = 0.003$$" with lines lc rgb "#0000FF"'\
',"$<" index 21 using 1:2 title "$$T/T_F = 0.060$$" with lines lc rgb "#006400"'\
',"$<" index 31 using 1:2 title "$$T/T_F = 0.090$$" with lines lc rgb "#32CD32"'\
',"$<" index 41 using 1:2 title "$$T/T_F = 0.130$$" with lines lc rgb "#FF8C00"'\
',"$<" index 51 using 1:2 title "$$T/T_F = 0.160$$" with lines lc rgb "#FF0000"'\
| gnuplot

Tdepend.tex: Tdepend Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$T/T_F$$"; 		   '\
'set ylabel "$$\max_k[\Delta_k]/\epsilon_F$$"; 	   '\
'plot [0:0.2][0:0.3] 		   '\
' "$<" index 0 using 1:3 notitle with lines lc rgb "black" '\
| gnuplot


kdepend Tdepend: Delta
	./$< 1> kdepend 2> Tdepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex Tdepend kdepend