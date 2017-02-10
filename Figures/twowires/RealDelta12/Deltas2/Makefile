CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: E0ddepend kdepend E0ddepend.tex kdepend.tex

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
' "$<" index 1     using 1:2 title "$$Interwire pairing real$$" with lines lc rgb "black"'\
| gnuplot

kdepend.tex: kdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$\Delta_k/\epsilon_{F,0}$$";'\
'set key top right;					'\
'plot  [-10:10][-0.6:0.6]			'\
' "$<" index 12  using 1:2 title "$$k_Fd = 0.720$$" with lines lc rgb "#e00000"'\
',"$<" index 12  using 1:3 notitle with lines lc rgb "#e00000"'\
',"$<" index 13  using 1:2 title "$$k_Fd = 0.735$$" with lines lc rgb "#e05900"'\
',"$<" index 13  using 1:3 notitle with lines lc rgb "#e05900"'\
',"$<" index 14  using 1:2 title "$$k_Fd = 0.750$$" with lines lc rgb "#167000"'\
',"$<" index 14  using 1:3 notitle with lines lc rgb "#167000"'\
',"$<" index 15  using 1:2 title "$$k_Fd = 0.765$$" with lines lc rgb "#002daa"'\
',"$<" index 15  using 1:3 notitle with lines lc rgb "#002daa"'\
',"$<" index 16  using 1:2 title "$$k_Fd = 0.775$$" with lines lc rgb "#9100aa"'\
',"$<" index 16  using 1:3 notitle with lines lc rgb "#9100aa"'\
| gnuplot


E0ddepend kdepend: Delta
	./$< 1> E0ddepend 2> kdepend

Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex *.gpi E0ddepend kdepend *.log