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
'set xlabel "$$k/kF$$"; 		   '\
'set ylabel "$$Delta/eF$$"; 	   '\
'plot [-10:10][-0.02:0.02] 		   '\
' "$<" index 0  using 1:2 notitle with lines '\
',"$<" index 40 using 1:2 notitle with lines '\
',"$<" index 60 using 1:2 notitle with lines '\
',"$<" index 80 using 1:2 notitle with lines '\
',"$<" index 90 using 1:2 notitle with lines '\
| gnuplot

Tdepend.tex: Tdepend Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$T/TF$$"; 		   '\
'set ylabel "$$Deltamax/eF$$"; 	   '\
'plot [0:0.01][0:0.02] 		   '\
' "$<" index 0 using 1:3 notitle with lines lc rgb "black" '\
| gnuplot


kdepend Tdepend: Delta
	./$< 1> kdepend 2> Tdepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf Tdepend kdepend