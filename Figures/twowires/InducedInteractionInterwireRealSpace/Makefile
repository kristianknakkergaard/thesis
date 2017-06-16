CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: InducedInteraction.tex retneg.tex data log InducedInteraction

InducedInteraction.tex: data makefile
	echo \
'set terminal epslatex color solid ; 								'\
'set output "$@"; 			  									 	'\
'set grid lw 2;       			   									'\
'set tics out;       			   									'\
'set xlabel "$$k_Fx$$"; 		   									'\
'set ylabel "$$V^{12}_{\\text{ind}}(x,0)/\\epsilon_{F,0}$$"; 	   	'\
'set key bottom left;												'\
'set key width -0.8;												'\
'set arrow from 0,0 to 0,-3 nohead lc rgb "black"; 					'\
'plot [-5.0:5.0][-3.0:0] 		   									'\
'0 notitle with lines lc rgb "black" 										'\
',"$<" index 5 using 1:2 title "$$k_Fd = 2.0$$" with lines lw 2 lc rgb "#002daa" '\
',"$<" index 4 using 1:2 title "$$k_Fd = 1.5$$" with lines lw 2 lc rgb "#167000" '\
',"$<" index 3 using 1:2 title "$$k_Fd = 1.0$$" with lines lw 2 lc rgb "#e05900" '\
',"$<" index 2 using 1:2 title "$$k_Fd = 0.5$$" with lines lw 2 lc rgb "#e00000" '\
',"$<" index 1 using 1:2 title "$$k_Fd = 0.0$$" with lines lw 2 lc rgb "black" '\
| gnuplot


retneg.tex: data makefile
	echo \
'set terminal epslatex color solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$n_F^{1/3}/n_B$$"; 		   '\
'set ylabel "$$retneg$$"; 	   '\
'plot [0:0.8][0:100.0] 		   '\
'0 notitle with lines lc rgb "black" '\
', 10.0 notitle with lines lc rgb "red" '\
',"$<" index 1 using 1:3 notitle with lines lc rgb "black" '\
| gnuplot


data log: InducedInteraction
	./$< 1> data 2> log  

InducedInteraction: InducedInteraction.o

clean:
	$(RM) log InducedInteraction InducedInteraction.o *.pdf *.eps *.tex data *.gpi