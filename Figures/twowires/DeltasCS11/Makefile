CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: ddepend.tex ddepend CS11depend.tex

CS11depend.tex: ddepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k_Fd$$"; 		    '\
'set ylabel "$$CS_{1,1}$$";         '\
'set key top right Left;				    '\
'set key width -6.5;				'\
'plot  [0.74:0.775][0:1.2]			'\
' "$<" index 0  using 1:6 title "$$Interwire pairing imaginary$$" with lines lc rgb "red"'\
',"$<" index 1  using 1:6 notitle with lines lc rgb "red"'\
',"$<" index 2  using 1:2 title "$$Interwire pairing real$$" with lines lc rgb "blue"'\
| gnuplot

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex ddepend