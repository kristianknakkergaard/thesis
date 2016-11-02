CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: findTC.tex findTC

findTC.tex: findTC Makefile
	echo \
'set terminal epslatex solid; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$T$$"; 		   '\
'set ylabel "$$Largest eigenvalue$$";   	   '\
'plot [0:][0:] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 1  using 1:2 notitle with lines lc rgb "blue"'\
',"$<" index 1  using 1:3 notitle with lines lc rgb "red"'\
| gnuplot

findTC: findTCcode
	./$< > $@

findTCcode: findTCcode.o

clean:
	$(RM) findTCcode findTCcode.o *.gpi *.pdf *.eps *depend.tex findTC