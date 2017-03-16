CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: xdepend11.tex xdepend12.tex xdepend 

xdepend11.tex: xdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k_Fx$$"; 		    '\
'set ylabel "$$Pair wave functions$$";'\
'set key bottom right;					'\
'plot [-15:15][-0.18:0.18]			'\
' "$<" index 3    using 1:2 title "$$k_Fd = 0.720$$" with lines lc rgb "#e00000"'\
',"$<" index 6    using 1:2 title "$$k_Fd = 0.735$$" with lines lc rgb "#e05900"'\
',"$<" index 9    using 1:2 title "$$k_Fd = 0.750$$" with lines lc rgb "#167000"'\
',"$<" index 12   using 1:2 title "$$k_Fd = 0.765$$" with lines lc rgb "#002daa"'\
',"$<" index 14   using 1:2 title "$$k_Fd = 0.775$$" with lines lc rgb "#9100aa"'\
| gnuplot

xdepend12.tex: xdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k_Fx$$"; 		    '\
'set ylabel "$$Pair wave functions$$";'\
'set key bottom right;					'\
'plot [-15:15][-0.18:0.18]			'\
' "$<" index 3    using 1:3 title "$$k_Fd = 0.720$$" with lines lc rgb "#e00000"'\
',"$<" index 6    using 1:3 title "$$k_Fd = 0.735$$" with lines lc rgb "#e05900"'\
',"$<" index 9    using 1:3 title "$$k_Fd = 0.750$$" with lines lc rgb "#167000"'\
',"$<" index 12   using 1:3 title "$$k_Fd = 0.765$$" with lines lc rgb "#002daa"'\
',"$<" index 14   using 1:3 title "$$k_Fd = 0.775$$" with lines lc rgb "#9100aa"'\
| gnuplot



xdepend: Delta
	./$< 1> xdepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex xdepend *.gpi