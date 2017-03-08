CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: corfunc.tex corfunc

corfunc.tex: corfunc Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k_Fx$$"; 		   '\
'set ylabel "$$ \braket{\psi_F(x)\psi_F(0)}/k_F $$"; 	   '\
'plot [-50:50][-0.1:0.1] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 2 using 1:2 title "$$T/T_F = 0.000$$" with lines lc rgb "#e00000" '\
',"$<" index 3 using 1:2 title "$$T/T_F = 0.075$$" with lines lc rgb "#e05900" '\
',"$<" index 4 using 1:2 title "$$T/T_F = 0.100$$" with lines lc rgb "#167000" '\
',"$<" index 5 using 1:2 title "$$T/T_F = 0.120$$" with lines lc rgb "#002daa" '\
',"$<" index 6 using 1:2 title "$$T/T_F = 0.127$$" with lines lc rgb "#9100aa" '\
| gnuplot


corfunc: Delta
	./$< 2> $@


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps corfunc.tex corfunc *.gpi