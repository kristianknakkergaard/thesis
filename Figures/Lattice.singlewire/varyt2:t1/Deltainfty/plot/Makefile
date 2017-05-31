CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: winding.tex

winding.tex: E0diffcorrected.txt Makefile
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

clean:
	$(RM) *.tex