CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: E0diff.txt E0diff.tex

E0diff.tex: E0diff.txt Makefile
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
' "$<" matrix with image notitle  	'\
| gnuplot

clean:
	$(RM) *.pdf *.eps *depend.tex xifilldepend *.gpi







