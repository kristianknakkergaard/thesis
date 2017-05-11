CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: xifilldepend xifilldepend.tex Delta

xifilldepend.tex: xifilldepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set palette defined (0 "white", 1 "red");							'\
'set xlabel "$$n$$"; 		    	'\
'set ylabel "$$\xi / d$$";			'\
'set xrange  [0:1]; 				'\
'set yrange  [1:11]; 				'\
'unset key;							'\
'unset colorbox;					'\
'plot 								'\
' "$<"  using 1:2:3 with image pixels notitle '\
| gnuplot

xifilldepend: Delta
	./$< 1> xifilldepend

Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex xifilldepend *.gpi