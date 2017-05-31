CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: E0diff.tex xifilldepend winding1k.tex winding2k.tex Delta

winding1k.tex: xifilldepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set palette defined (0 "white", 1 "red" );	'\
'set xlabel "$$n$$"; 		    	'\
'set ylabel "$$\xi / d$$";			'\
'unset key;							'\
'unset colorbox;					'\
'set xrange  [0:1]; 				'\
'set yrange  [1:7]; 				'\
'plot 								'\
' "$<"  using 1:2:3 with image pixels notitle '\
| gnuplot

winding2k.tex: xifilldepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set palette defined (0 "white", 1 "red", 2 "blue" ); '\
'set xlabel "$$n$$"; 		    	'\
'set ylabel "$$\xi / d$$";			'\
'unset key;							'\
'unset colorbox;					'\
'set xrange  [0:1]; 				'\
'set yrange  [1:7]; 				'\
'plot 								'\
' "$<"  using 1:2:4 with image pixels notitle '\
| gnuplot

E0diff.tex: xifilldepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set autoscale xfix; 				'\
'set autoscale yfix; 				'\
'set palette defined ( 0 "white", 1 "blue");	'\
'set xlabel "$$n$$"; 		    	'\
'set ylabel "$$\xi / d$$";			'\
'unset key;							'\
'set xrange  [0:1]; 				'\
'set yrange  [1:7]; 				'\
'set colorbox;						'\
'plot 								'\
' "$<" using 1:2:($$6 - $$5) with image notitle  	'\
| gnuplot

xifilldepend: Delta
	./$< 1> xifilldepend

Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *.tex xifilldepend