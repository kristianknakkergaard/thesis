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
'plot [-10:10][-0.6:0.6]			'\
' "$<" index 3    using 1:2 title "$$k_Fd = 0.575$$" with lines lc rgb "#e00000"'\
',"$<" index 3    using 1:3 notitle with lines lc rgb "#e00000"'\
',"$<" index 5    using 1:2 title "$$k_Fd = 0.585$$" with lines lc rgb "#e05900"'\
',"$<" index 5    using 1:3 notitle with lines lc rgb "#e05900"'\
',"$<" index 7    using 1:2 title "$$k_Fd = 0.595$$" with lines lc rgb "#167000"'\
',"$<" index 7    using 1:3 notitle with lines lc rgb "#167000"'\
',"$<" index 10   using 1:2 title "$$k_Fd = 0.610$$" with lines lc rgb "#002daa"'\
',"$<" index 10   using 1:3 notitle with lines lc rgb "#002daa"'\
',"$<" index 12   using 1:2 title "$$k_Fd = 0.620$$" with lines lc rgb "#9100aa"'\
',"$<" index 12   using 1:3 notitle with lines lc rgb "#9100aa"'\
| gnuplot



kdepend: Delta
	./$< 1> kdepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex kdepend *.gpi