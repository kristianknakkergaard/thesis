CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend.tex Tdepend.tex chempot.tex kdepend Tdepend

kdepend.tex: kdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$\Delta_k/\epsilon_{F,0}$$";'\
'set key top left;					'\
'plot [-10:10][-0.35:0.35]			'\
' "$<" index 0  using 1:2 title  "$$T/T_F = 0.000$$" with lines lc rgb "#e00000"'\
',"$<" index 100 using 1:2 title "$$T/T_F = 0.100$$" with lines lc rgb "#e05900"'\
',"$<" index 130 using 1:2 title "$$T/T_F = 0.130$$" with lines lc rgb "#167000"'\
',"$<" index 145 using 1:2 title "$$T/T_F = 0.145$$" with lines lc rgb "#002daa"'\
| gnuplot

Tdepend.tex: Tdepend Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set key top right;       	   '\
'set xlabel "$$T/T_F$$"; 		   '\
'set ylabel "$$\max_k[\Delta_k]/\epsilon_{F,0}$$"; 	   '\
'plot [0:0.2][0:0.35] 		   '\
' "$<" index 0 using 1:3 title "Numerical calculation" with lines lc rgb "black" '\
',"$<" index 1 using 1:2 title "Asymptotic form" with lines lc rgb "red" '\
| gnuplot

chempot.tex: Tdepend Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set key top left;       	   '\
'set xlabel "$$T/T_F$$"; 		   '\
'set ylabel "$$\mu/\epsilon_{F,0}$$"; 	   '\
'plot [0:0.2][0.94:1.04] 		   '\
' "$<" index 0 using 1:5 title "Free gas" with lines lc rgb "red" '\
',"$<" index 0 using 1:4 title "Interacting gas" with lines lc rgb "black" '\
| gnuplot


kdepend Tdepend: Delta
	./$< 1> kdepend 2> Tdepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *.tex *.gpi Tdepend kdepend