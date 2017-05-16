CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: ddepend.tex ddepend CS11depend.tex

ddepend.tex: ddepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k_Fd$$"; 		    '\
'set ylabel "$$Delta_k/epsilon_{F,0}$$";'\
'set key top left;					'\
'set key width -6.5;				'\
'plot [0.71:0.78][0:0.6]			'\
' "$<" index 1  using 1:2 title "$$Intrawire pairing$$" with lines lc rgb "red"'\
',"$<" index 1  using 1:3 title "$$Interwire pairing$$" with lines lc rgb "blue"'\
',"$<" index 2  using 1:2 notitle with lines lc rgb "red"'\
',"$<" index 2  using 1:3 notitle with lines lc rgb "blue"'\
| gnuplot

CS11depend.tex: ddepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k_Fd$$"; 		    '\
'set ylabel "$$CS_{1,1}$$";         '\
'set key top right;				    '\
'set key width -6.5;				'\
'plot  [0.71:0.78][0:1.2]			'\
' "$<" index 1  using 1:6 title "$$Interwire pairing imaginary$$" with lines lc rgb "red"'\
',"$<" index 2  using 1:6 notitle with lines lc rgb "red"'\
',"$<" index 3  using 1:2 title "$$Interwire pairing real$$" with lines lc rgb "blue"'\
| gnuplot


ddepend: Delta
	./$< > ddepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex ddepend