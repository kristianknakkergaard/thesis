CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: ddepend.tex mudepend.tex Edepend.tex ddepend

ddepend.tex: ddepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k_Fd$$"; 		    '\
'set ylabel "$$Delta_k/epsilon_{F,0}$$";'\
'set key top right Left;			'\
'set key width -4.5;				'\
'plot  [0.74:0.775][0:0.4]			'\
' "$<" index 1  using 1:2 title "$$Intrawire pairing 1$$" with lines lc rgb "#4a0000"'\
',"$<" index 1  using 1:3 title "$$Interwire pairing 1$$" with lines dt 4 lc rgb "#ff0000"'\
',"$<" index 2  using 1:2 notitle with lines lc rgb "#4a0000"'\
',"$<" index 2  using 1:3 notitle with lines dt 4 lc rgb "#ff0000"'\
| gnuplot

Edepend.tex: ddepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k_Fd$$"; 		    '\
'set ylabel "$$E_0 / 2\epsilon_{F,0}N_F$$";'\
'set key top left;					'\
'plot [0.74:0.775][0.65:0.655]			'\
' "$<" index 1  using 1:($$4 + 2 * $$5) notitle with lines lc rgb "red"'\
',"$<" index 2  using 1:($$4 + 2 * $$5) notitle with lines lc rgb "red"'\
| gnuplot

mudepend.tex: ddepend Makefile
	echo \
'set terminal epslatex solid ; '\
'set output "$@"; 			   '\
'set grid;       			   '\
'set tics out;       			   '\
'set xlabel "$$k_Fd$$"; 		   '\
'set ylabel "$$\mu(0)/\epsilon_{F,0}"; 	   '\
'plot  [0.74:0.775][0.95:1.1] 		   '\
' "$<" index 1 using  1:5 notitle with lines lc rgb "black" '\
',"$<" index 2 using  1:5 notitle with lines lc rgb "black" '\
',"$<" index 1 using  1:(1.0) notitle with lines lc rgb "red" '\
',"$<" index 2 using  1:(1.0) notitle with lines lc rgb "red" '\
| gnuplot

ddepend: Delta
	./$< 1> ddepend

Delta: Delta.o

clean:
	$(RM) *.pdf *.eps *depend.tex Delta *.o