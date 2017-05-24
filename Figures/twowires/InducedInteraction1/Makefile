CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: InducedInteraction.tex data log InducedInteraction

InducedInteraction.tex: data makefile
	echo \
'set terminal epslatex color solid ; 	'\
'set output "$@"; 			   			'\
'set grid;       			   			'\
'set xlabel "$$n_F / n_B^{1/3}$$"; 		'\
'set ylabel "$$k_Fd_c$$"; 	   			'\
'plot [0:1.5][0:1.2] 				    '\
' "$<" index 0 using 1:2 with lines notitle '\
| gnuplot

data log: InducedInteraction
	./$< 1> data 2> log  

InducedInteraction: InducedInteraction.o

clean:
	$(RM) log InducedInteraction InducedInteraction.o *.pdf *.eps *.tex data *.gpi