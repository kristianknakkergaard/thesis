CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: E0depend.tex xifilldepend

E0depend.tex: xifilldepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set autoscale xfix; 				'\
'set autoscale yfix; 				'\
'set palette grey;	'\
'set xlabel "$$n$$"; 		    	'\
'set ylabel "$$\xi / d$$";			'\
'unset key;							'\
'set xrange  [0:1]; 				'\
'set yrange  [1:7]; 				'\
'set colorbox;						'\
'plot 								'\
' "$<" using 1:2:5 with image notitle  	'\
| gnuplot

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex xifilldepend *.gpi