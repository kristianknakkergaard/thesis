CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend.tex nBdepend.tex kdepend parameters Delta

kdepend.tex: kdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$Delta_k/epsilon_{F,0}$$";'\
'set key bottom right;					'\
'set arrow from 0,-3 to 0,3 nohead lc rgb "black"; '\
'set arrow from -3.2,0 to 3.2,0 nohead lc rgb "black"; '\
'plot [-3.14:3.14][-3:3]			'\
' "$<" index 45  using 1:2 title "$$only \sin(k)$$" with lines lc rgb "#e00000"'\
',"$<" index 45  using 1:3 title "$$only \sin(2k)$$" with lines lc rgb "blue"'\
',"$<" index 45  using 1:4 title "$$\varepsilon_k$$" with lines lc rgb "black"'\
| gnuplot

nBdepend.tex: parameters Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$n_B/n_F^3$$"; 		    '\
'set ylabel "$$Ground state energy$$";'\
'set key bottom right;					'\
'set arrow from 0,-3 to 0,0 nohead lc rgb "black"; '\
'set arrow from 0,0 to 2,0 nohead lc rgb "black"; '\
'plot  [0:0.0001][-6.15:-5.6]			'\
' "$<" index 0  using 1:2 title "$$only \sin(k)$$" with lines lc rgb "#e00000"'\
',"$<" index 0  using 1:3 title "$$only \sin(2k)$$" with lines lc rgb "blue"'\
| gnuplot

kdepend parameters: Delta
	./$< 1> kdepend 2> parameters


Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex parameters kdepend *.gpi