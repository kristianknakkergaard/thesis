CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: ndepend ndepend.tex Delta

ndepend.tex: ndepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$n$$"; 		  		'\
'set ylabel "$$\mu$$";				'\
'unset key;							'\
'set arrow from 0,0 to 1,0 nohead lc rgb "black"; '\
'plot [0:1][-3.0:3.0]			'\
' "$<" index 0  using 1:2 notitle with lines lc rgb "red"'\
',"$<" index 0  using 1:3 notitle with lines lc rgb "blue"'\
| gnuplot

ndepend: Delta
	./$< 1> ndepend

Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex ndepend *.gpi