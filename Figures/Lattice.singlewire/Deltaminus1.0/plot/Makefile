CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: winding.tex E0diff.tex

winding.tex: Deltamaxcorrected.txt Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set palette defined (0 "white", 1 "red", 2 "blue" );	'\
'set xlabel "$$n$$"; 		    	'\
'set ylabel "$$\xi / d$$";			'\
'unset key;							'\
'unset colorbox;					'\
'set xrange  [0:1]; 				'\
'set yrange  [1:7]; 				'\
'plot 								'\
' "$<"  using 1:2:3 with image pixels notitle '\
| gnuplot

E0diff.tex: Deltamaxcorrected.txt Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set autoscale xfix; 				'\
'set autoscale yfix; 				'\
'set palette defined (0 "white", 1 "blue");	'\
'set xlabel "$$n$$"; 		    	'\
'set ylabel "$$\xi / d$$";			'\
'unset key;							'\
'set xrange  [0:1]; 				'\
'set yrange  [1:7]; 				'\
'set colorbox;						'\
'set cbrange [0:0.03];				'\
'plot 								'\
' "$<" using 1:2:4 with image notitle  	'\
| gnuplot

clean:
	$(RM) *.tex