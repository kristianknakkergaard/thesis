CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: Deltamax.tex xifilldepend xifilldepend.tex Delta

xifilldepend.tex: xifilldepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set palette defined (0 "white", 1 "red", 2 "blue" );							'\
'set xlabel "$$n$$"; 		    	'\
'set ylabel "$$\xi / d$$";			'\
'unset key;							'\
'unset colorbox;					'\
'set xrange  [0:1]; 				'\
'set yrange  [1:11]; 				'\
'plot 								'\
' "$<"  using 1:2:3 with image pixels notitle '\
| gnuplot

Deltamax.tex: xifilldepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set autoscale xfix; 				'\
'set autoscale yfix; 				'\
'set palette defined ( 0 "green", 1 "blue", 2 "red", 3 "orange" );	'\
'set xlabel "$$n$$"; 		    	'\
'set ylabel "$$\xi / d$$";			'\
'unset key;							'\
'set colorbox;						'\
'plot 								'\
' "$<" using 1:2:5 with image notitle  	'\
| gnuplot

xifilldepend: Delta
	./$< 1> xifilldepend

Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex xifilldepend *.gpi