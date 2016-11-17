CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: TCeigen.tex TCeigen.txt

TCeigen.tex: TCeigen.txt Makefile
	echo \
'set terminal epslatex solid; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$T/T_F$$"; 		   '\
'set ylabel "$$Eigenvalue$$";   	   '\
'plot [0:0.2][0:2] 		   '\
'0 notitle with lines lc rgb "black" '\
', 1 notitle with lines lc rgb "red" '\
',"$<" index 0  using 1:2 notitle with lines lc rgb "black"'\
',"$<" index 0  using 1:3 notitle with lines lc rgb "black"'\
',"$<" index 0  using 1:4 notitle with lines lc rgb "black"'\
',"$<" index 0  using 1:5 notitle with lines lc rgb "black"'\
',"$<" index 0  using 1:6 notitle with lines lc rgb "black"'\
| gnuplot

clean:
	$(RM) *.gpi *.pdf *.eps *TCeigen.tex