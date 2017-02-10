CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: E0ddepend E0ddepend.tex

E0ddepend.tex: E0ddepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k_Fd$$"; 		    '\
'set ylabel "$$(E_0 + \mu N_F)/(\varepsilon_{F,0} N_F)$$";'\
'set key top right;					'\
'plot [0.74:0.75][0.61:0.69]			'\
' "$<" index 1     using 1 : (2.0 / 3.0) title "$$Free gas$$" with lines lw 2 lc rgb "black"  '\
',"$<" index 1     using 1:(0.628777) title "$$Only intrawire pairing$$" with lines dt 2 lw 2 lc rgb "black"'\
',"$<" index 1     using 1:3 title "$$Only interwire pairing$$" with lines dt 4 lw 2 lc rgb "black" '\
',"$<" index 1     using 1:2 title "$$Interwire pairing imaginary$$" with lines lc rgb "blue"'\
| gnuplot

E0ddepend: Delta
	./$< 1> E0ddepend

Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex *.gpi E0ddepend