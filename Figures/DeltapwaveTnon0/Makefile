CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend.tex Tdepend.tex kdepend Tdepend

kdepend.tex: kdepend Makefile
	echo \
'set terminal epslatex; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k/k_F$$"; 		   '\
'set ylabel "$$Delta_k/epsilon_F$$";   	   '\
'set key top left;		   	   '\
'plot [-20:20][-0.25:0.25] 		   '\
' "$<" index 0  using 1:2 title "T/T_F = 0.00" with lines lc rgb "#e00000"'\
',"$<" index 40 using 1:2 title "T/T_F = 0.08" with lines lc rgb "#e05900"'\
',"$<" index 50 using 1:2 title "T/T_F = 0.10" with lines lc rgb "#167000"'\
',"$<" index 60 using 1:2 title "T/T_F = 0.12" with lines lc rgb "#002daa"'\
',"$<" index 70 using 1:2 title "T/T_F = 0.14" with lines lc rgb "#9100aa"'\
| gnuplot

Tdepend.tex: Tdepend Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$T/T_F$$"; 		   '\
'set ylabel "$$\max_k[\Delta_k]/\epsilon_F$$"; 	   '\
'plot [0:0.2][0:0.25] 		   '\
' "$<" index 0 using 1:3 notitle with lines lc rgb "black" '\
| gnuplot


kdepend Tdepend: Delta
	./$< 1> kdepend 2> Tdepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex Tdepend kdepend