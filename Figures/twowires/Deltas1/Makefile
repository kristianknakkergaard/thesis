CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend.tex kdepend 

kdepend.tex: kdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$Delta_k/epsilon_{F,0}$$";'\
'set key bottom right;					'\
'set key width -1;						'\
'set arrow from 0,-0.6 to 0,0.6 nohead lc rgb "black"; 	'\
'set arrow from -10,0 to 10,0 nohead lc rgb "black"; 	'\
'plot [-10:10][-0.6:0.6]			'\
' "$<" index 3    using 1:2 title "$$k_Fd = 0.720$$" with lines lc rgb "#e00000"'\
',"$<" index 3    using 1:3 notitle with lines lc rgb "#e00000"'\
',"$<" index 6    using 1:2 title "$$k_Fd = 0.735$$" with lines lc rgb "#e05900"'\
',"$<" index 6    using 1:3 notitle with lines lc rgb "#e05900"'\
',"$<" index 9    using 1:2 title "$$k_Fd = 0.750$$" with lines lc rgb "#167000"'\
',"$<" index 9    using 1:3 notitle with lines lc rgb "#167000"'\
',"$<" index 12   using 1:2 title "$$k_Fd = 0.765$$" with lines lc rgb "#002daa"'\
',"$<" index 12   using 1:3 notitle with lines lc rgb "#002daa"'\
',"$<" index 14   using 1:2 title "$$k_Fd = 0.775$$" with lines lc rgb "#9100aa"'\
',"$<" index 14   using 1:3 notitle with lines lc rgb "#9100aa"'\
| gnuplot


kdepend: Delta
	./$< 1> kdepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex kdepend *.gpi