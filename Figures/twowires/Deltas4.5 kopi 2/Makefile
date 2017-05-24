CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: Tdepend.tex mudepend.tex Tdepend

Tdepend.tex: Tdepend Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$T/T_F$$"; 		   '\
'set ylabel "$$\max_k[\Delta_k]/\epsilon_{F,0}$$"; 	   '\
'plot [0:0.3][0:0.5] 		   '\
' "$<" index 1  using 1:2 notitle with lines lc rgb "red" '\
',"$<" index 2  using 1:(1.3 * $$2) notitle with lines dt 2 lc rgb "red" '\
',"$<" index 1  using 1:3 notitle with lines lc rgb "blue" '\
',"$<" index 3  using 1:2 notitle with lines dt 2 lc rgb "blue" '\
| gnuplot

mudepend.tex: Tdepend Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$T/T_F$$"; 		   '\
'set ylabel "$$\mu(T)/\epsilon_{F,0}"; 	   '\
'plot [0:0.5][0.95:1.1] 		   '\
' "$<"  index 1 using 1:4 notitle with lines lc rgb "black" '\
',"$<"  index 1 using 1:(1.0 + 3.14159265359 * 3.14159265359 / 12.0 * $$1 * $$1) notitle with lines lc rgb "red" '\
| gnuplot


Tdepend: Delta
	./$< 1> Tdepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex Tdepend kdepend *.gpi