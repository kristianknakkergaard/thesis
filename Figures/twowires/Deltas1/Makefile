CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend.tex Tdepend.tex mudepend.tex kdepend Tdepend

kdepend.tex: kdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$Delta_k/epsilon_{F,0}$$";'\
'set key top left;					'\
'plot [-10:10][-0.3:0.3]			'\
' "$<" index 0     using 1:2 title "$$\Delta^{11}_k$$" with lines lc rgb "red"'\
',"$<" index 0     using 1:3 title "$$\Delta^{12}_k$$" with lines lc rgb "blue"'\
#',"$<" index 750   using 1:2 title "$$T/T_F = 0.075$$" with lines lc rgb "#e05900"'\
#',"$<" index 1000  using 1:2 title "$$T/T_F = 0.100$$" with lines lc rgb "#167000"'\
#',"$<" index 1200  using 1:2 title "$$T/T_F = 0.120$$" with lines lc rgb "#002daa"'\
#',"$<" index 1270  using 1:2 title "$$T/T_F = 0.127$$" with lines lc rgb "#9100aa"'\
| gnuplot

Tdepend.tex: Tdepend Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$T/T_F$$"; 		   '\
'set ylabel "$$\max_k[\Delta_k]/\epsilon_{F,0}$$"; 	   '\
'plot [0:0.15][0:0.3] 		   '\
' "$<" index 1 using 1:3 notitle with lines lc rgb "black" '\
',"$<" index 2 using 1:2 notitle with lines lc rgb "red" '\
| gnuplot


mudepend.tex: Tdepend Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$T/T_F$$"; 		   '\
'set ylabel "$$\mu(T)/\epsilon_{F,0}"; 	   '\
'plot [0:0.15][0.95:1.1] 		   '\
' "$<" index 1 using 1:4 notitle with lines lc rgb "black" '\
',"$<" index 1 using 1:(1.0 + 3.14159265359 * 3.14159265359 / 12.0 * $$1 * $$1) notitle with lines lc rgb "red" '\
| gnuplot


kdepend Tdepend: Delta
	./$< 1> kdepend 2> Tdepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex Tdepend kdepend *.gpi