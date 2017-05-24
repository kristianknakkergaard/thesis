CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: Occupancyplot.tex energyplot.tex data Occupancy

Occupancyplot.tex: data makefile
	echo \
'set terminal epslatex color solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k/k_F$$"; 		   '\
'set ylabel "$$\\braket{c_{1,k}^\\dagger c_{1,k}}$$"; 	   '\
'set key top right Left;					'\
'set key width -2;'\
'plot [0:2][0:1.1] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 0 using 1:2 title "Free gas" with lines lw 2 lc rgb "black" '\
',"$<" index 1 using 1:2 title "$$k_Fd = 0.70$$" with lines lw 2 lc rgb "#e00000" '\
',"$<" index 2 using 1:2 title "$$k_Fd = 0.75$$" with lines lw 2 lc rgb "#167000" '\
',"$<" index 3 using 1:2 title "$$k_Fd = 0.80$$" with lines lw 2 lc rgb "#9100aa" '\
| gnuplot


energyplot.tex: data makefile
	echo \
'set terminal epslatex color solid ; 		'\
'set output "$@"; 			   				'\
'set grid;       			   				'\
'set tics out;       			   			'\
'set xlabel "$$k/k_F$$"; 		   			'\
'set ylabel "$$E_{F,k}/\\epsilon_{F,0}$$"; 	'\
'set key top right Left;					'\
'set key width -2;							'\
'plot [0:2.0][0:5.0] 		   				'\
'0 notitle with lines lc rgb "black" '\
',"$<" index 0 using 1:5 title "Free gas" with lines lw 2 lc rgb "black" '\
',"$<" index 1 using 1:4 title "$$k_Fd = 0.70$$" with lines lw 2 lc rgb "#e00000" '\
',"$<" index 2 using 1:4 title "$$k_Fd = 0.75$$" with lines lw 2 lc rgb "#167000" '\
',"$<" index 3 using 1:4 title "$$k_Fd = 0.80$$" with lines lw 2 lc rgb "#9100aa" '\
| gnuplot

data: Occupancy
	./$< > data  

Occupancy: Occupancy.o

clean:
	$(RM) log Occupancy Occupancy.o *.pdf *.eps *plot.tex data *.gpi