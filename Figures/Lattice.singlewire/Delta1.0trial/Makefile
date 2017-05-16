CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: xifilldepend.tex xifilldepend

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
'plot  [0:1][1:7]								'\
' "$<"  using 1:2:3 with image pixels notitle '\
| gnuplot

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex xifilldepend *.gpi







