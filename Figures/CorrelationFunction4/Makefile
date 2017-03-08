CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: corfunc.tex aveordparam.tex corfunc params

corfunc.tex: corfunc Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k/k_F$$"; 		   '\
'set ylabel "$$Correlation function$$"; 	   '\
'plot [0:2][-0.5:0.5] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 2 using 1:2 title "$$T/T_F = 0.000$$" with lines lc rgb "#e00000" '\
',"$<" index 3 using 1:2 title "$$T/T_F = 0.075$$" with lines lc rgb "#e05900" '\
',"$<" index 4 using 1:2 title "$$T/T_F = 0.100$$" with lines lc rgb "#167000" '\
',"$<" index 5 using 1:2 title "$$T/T_F = 0.120$$" with lines lc rgb "#002daa" '\
',"$<" index 6 using 1:2 title "$$T/T_F = 0.127$$" with lines lc rgb "#9100aa" '\
| gnuplot

aveordparam.tex: corfunc Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k / k_F $$"; 		   '\
'set ylabel "$$ Mean field squared $$"; 	   '\
'plot [0:2][-0.5:0.5] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 2 using 1:3 title "$$T/T_F = 0.000$$" with lines lc rgb "#e00000" '\
',"$<" index 3 using 1:3 title "$$T/T_F = 0.075$$" with lines lc rgb "#e05900" '\
',"$<" index 4 using 1:3 title "$$T/T_F = 0.100$$" with lines lc rgb "#167000" '\
',"$<" index 5 using 1:3 title "$$T/T_F = 0.120$$" with lines lc rgb "#002daa" '\
',"$<" index 6 using 1:3 title "$$T/T_F = 0.127$$" with lines lc rgb "#9100aa" '\
| gnuplot


corfunc params: Delta
	./$< 1> corfunc 2> params


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps corfunc.tex aveordparam.tex corfunc *.gpi