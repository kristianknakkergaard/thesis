CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: InducedInteraction.tex retneg.tex data log InducedInteraction

InducedInteraction.tex: data makefile
	echo \
'set terminal epslatex color solid ; 						'\
'set output "$@"; 			   								'\
'set grid;       			   								'\
'set tics out;       			   							'\
'set xlabel "$$q/k_F$$"; 		   							'\
'set ylabel "$$2m_F/k_F V_{FF,12}^\text{ind}(q,0)$$"; 	   	'\
'set key bottom left;										'\
'set key width -0.9;										'\
'set arrow from 0,0 to 0,-5 nohead lc rgb "black"; 			'\
'plot [-5.0:5.0][-5:0] 		   								'\
'0 notitle with lines lc rgb "black" 						'\
',"$<" index 4 using 1:2 title "$$k_Fd = 1.8$$" with lines lc rgb "#002daa" '\
',"$<" index 3 using 1:2 title "$$k_Fd = 1.4$$" with lines lc rgb "#167000" '\
',"$<" index 2 using 1:2 title "$$k_Fd = 1.0$$" with lines lc rgb "#e05900" '\
',"$<" index 1 using 1:2 title "$$k_Fd = 0.6$$" with lines lc rgb "#e00000" '\
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