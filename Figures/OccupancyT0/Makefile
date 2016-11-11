CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: Occupancyplot.tex energyplot.tex data log Occupancy

Occupancyplot.tex: data makefile
	echo \
'set terminal epslatex color solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k/k_F$$"; 		   '\
'set ylabel "$$\braket{f_k^\dagger f_k}$$"; 	   '\
'set key top right;					'\
'plot [0:2][0:1.1] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 0 using 1:2 title "Free gas" with lines lw 2 lc rgb "black" '\
',"$<" index 1 using 1:2 title "$$r_B = 0.005$$" with lines lw 2 lc rgb "#e00000" '\
',"$<" index 2 using 1:2 title "$$r_B = 0.006$$" with lines lw 2 lc rgb "#167000" '\
',"$<" index 3 using 1:2 title "$$r_B = 0.007$$" with lines lw 2 lc rgb "#002daa" '\
',"$<" index 4 using 1:2 title "$$r_B = 0.010$$" with lines lw 2 lc rgb "#9100aa" '\
| gnuplot


energyplot.tex: data makefile
	echo \
'set terminal epslatex color solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k/k_F$$"; 		   '\
'set ylabel "$$ E_{F,k}$$"; 	   '\
'set key top left;					'\
'plot [0:2.0][0:5.0] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 0 using 1:5 title "Free gas" with lines lw 2 lc rgb "black" '\
',"$<" index 1 using 1:5 title "$$r_B = 0.005$$" with lines lw 2 lc rgb "#e00000" '\
',"$<" index 2 using 1:5 title "$$r_B = 0.006$$" with lines lw 2 lc rgb "#167000" '\
',"$<" index 3 using 1:5 title "$$r_B = 0.007$$" with lines lw 2 lc rgb "#002daa" '\
',"$<" index 4 using 1:5 title "$$r_B = 0.010$$" with lines lw 2 lc rgb "#9100aa" '\
| gnuplot

data log: Occupancy
	./$< 1> data 2> log  

Occupancy: Occupancy.o

clean:
	$(RM) log Occupancy Occupancy.o *.pdf *.eps *plot.tex data *.gpi