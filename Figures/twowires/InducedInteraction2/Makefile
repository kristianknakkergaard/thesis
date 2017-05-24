CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: InducedInteraction.tex data log InducedInteraction

InducedInteraction.tex: data makefile
	echo \
'set terminal epslatex color solid ; 	'\
'set output "$@"; 			   			'\
'set grid;       			   			'\
'set tics out;       			   		'\
'set key top right Left; 				'\
'set key width -2; 						'\
'set key opaque; 						'\
'set xlabel "$$k / k_F$$"; 		'\
'set ylabel "$$2m_F/k_F W_{\text{ind}}(k, k_F)$$"; 	   '\
'plot [-10:10][-5:5] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 0   using 1:2 title "$$k_Fd = 0.720$$" with lines lc rgb "#e00000"'\
',"$<" index 0   using 1:3 notitle with lines lc rgb "#e00000"'\
',"$<" index 1   using 1:2 title "$$k_Fd = 0.735$$" with lines lc rgb "#e05900"'\
',"$<" index 1   using 1:3 notitle with lines lc rgb "#e05900"'\
',"$<" index 2   using 1:2 title "$$k_Fd = 0.750$$" with lines lc rgb "#167000"'\
',"$<" index 2   using 1:3 notitle with lines lc rgb "#167000"'\
',"$<" index 3   using 1:2 title "$$k_Fd = 0.765$$" with lines lc rgb "#002daa"'\
',"$<" index 3   using 1:3 notitle with lines lc rgb "#002daa"'\
',"$<" index 4   using 1:2 title "$$k_Fd = 0.775$$" with lines lc rgb "#9100aa"'\
',"$<" index 4   using 1:3 notitle with lines lc rgb "#9100aa"'\
| gnuplot

data log: InducedInteraction
	./$< 1> data 2> log  

InducedInteraction: InducedInteraction.o

clean:
	$(RM) log InducedInteraction InducedInteraction.o *.pdf *.eps *.tex data *.gpi