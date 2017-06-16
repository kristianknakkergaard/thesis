CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: xdepend11.tex xdepend12.tex xdepend 

xdepend11.tex: xdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid lw 2;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k_Fx$$"; 		    '\
'set ylabel "$$\braket{\psi_{1,F}(x)\psi_{1,F}(0)} / k_F$$";	'\
'set key bottom right;					'\
'set key width -1;						'\
'set arrow from 0,-0.12 to 0,0.12 nohead lc rgb "black"; 	'\
'set arrow from -15,0 to 15,0 nohead lc rgb "black"; 		'\
'plot [-15:15][-0.12:0.12]									'\
' "$<" index 4    using 1:2 title "$$k_Fd = 0.7475$$" with lines lw 2 lc rgb "#e00000"'\
',"$<" index 6    using 1:2 title "$$k_Fd = 0.7525$$" with lines lw 2 lc rgb "#e05900"'\
',"$<" index 8    using 1:2 title "$$k_Fd = 0.7575$$" with lines lw 2 lc rgb "#167000"'\
',"$<" index 10   using 1:2 title "$$k_Fd = 0.7625$$" with lines lw 2 lc rgb "#002daa"'\
',"$<" index 12   using 1:2 title "$$k_Fd = 0.7675$$" with lines lw 2 lc rgb "#9100aa"'\
| gnuplot

xdepend12.tex: xdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid lw 2;       			   	'\
'set tics out;       			    '\
'set xlabel "$$k_Fx$$"; 		    '\
'set ylabel "$$\braket{\psi_{1,F}(x)\psi_{2,F}(0)} / k_F$$";	'\
'set key bottom right;					'\
'set key width -1;						'\
'set arrow from 0,-0.12 to 0,0.12 nohead lc rgb "black"; '\
'set arrow from -15,0 to 15,0 nohead lc rgb "black"; '\
'plot [-15:15][-0.12:0.12]			'\
' "$<" index 4    using 1:3 title "$$k_Fd = 0.7475$$" with lines lw 2 lc rgb "#e00000"'\
',"$<" index 6    using 1:3 title "$$k_Fd = 0.7525$$" with lines lw 2 lc rgb "#e05900"'\
',"$<" index 8    using 1:3 title "$$k_Fd = 0.7575$$" with lines lw 2 lc rgb "#167000"'\
',"$<" index 10   using 1:3 title "$$k_Fd = 0.7625$$" with lines lw 2 lc rgb "#002daa"'\
',"$<" index 12   using 1:3 title "$$k_Fd = 0.7675$$" with lines lw 2 lc rgb "#9100aa"'\
| gnuplot



xdepend: Delta
	./$< 1> xdepend


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex xdepend *.gpi