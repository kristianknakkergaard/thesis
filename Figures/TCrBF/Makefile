CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: TCrBF.tex TCrBF.txt

TCrBF.tex: TCrBF.txt Makefile
	echo \
'set terminal epslatex solid; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$T/T_F$$"; 		   '\
'set ylabel "$$(n_Ba_{BF}^3)^{1/3}$$";   	   '\
'plot [0:0.1][0:0.13] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 1  using 1:2 notitle with lines lc rgb "black"'\
| gnuplot

clean:
	$(RM) *.gpi *.pdf *.eps *TCrBF.tex