CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: ddepend.tex ddepend

ddepend.tex: ddepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$Delta_k/epsilon_{F,0}$$";'\
'set key top left;					'\
'plot [0.59:0.615][0:0.5]			'\
' "$<" index 0    	using 1:2 title "$$\Delta^{11}_k$$" with lines lc rgb "red"'\
',"$<" index 0    	using 1:3 title "$$\Delta^{11}_k$$" with lines lc rgb "blue"'\
| gnuplot



ddepend: Delta
	./$< > ddepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex ddepend