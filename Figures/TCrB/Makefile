CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: TCrB.tex TCrBB.txt

TCrB.tex: TCrBB.txt Makefile
	echo \
'set terminal epslatex solid; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$T/T_F$$"; 		   '\
'set ylabel "$$(n_Ba_B^3)^{1/3}$$";   	   '\
'plot [0.005:0.03][0:0.25] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 1  using 1:2 notitle with lines lc rgb "black"'\
| gnuplot

clean:
	$(RM) *.gpi *.pdf *.eps *TCrB.tex