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
'plot [-20:20][0:5e-8] 		   '\
'0 notitle with lines lc rgb "black" '\
',"$<" index 2 using 1:2 notitle with lines lc rgb "#e00000" '\
| gnuplot


corfunc: Delta
	./$< 2> $@


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps corfunc.tex corfunc *.gpi