CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: nBdepend.tex nBdepend

nBdepend.tex: nBdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set termoption enhanced; 			'\
'set output "$@"; 			   		'\
'set grid front;       			   		'\
'set tics out;       			    '\
'set xlabel "$$n_F/n_B^{1/3}$$";    '\
'set ylabel "$$k_Fd_c$$";			'\
'set style fill pattern 3;			    				'\
'set label 1 at 0.22, 0.9 front;			                    '\
'set label 1 "$$Intrawire pairing only$$";			    '\
'set label 2 at 0.82, 0.5 front;			                    '\
'set label 2 "$$Interwire pairing only$$";			    '\
'plot [0:1.5][0:1.2]			    '\
' "$<" index 0  using 1:2 notitle with filledcurves x1 lw 2 lc rgb "grey"'\
',"$<" index 0  using 1:2 notitle with lines lw 2 lc rgb "black"'\
',"$<" index 0  using ($0.6447880954):2 notitle with lines dt 2 lw 2 lc rgb "black"'\
',"$<" index 5  using 1:2 notitle with lines dt 4 lw 2 lc rgb "black"'\
| gnuplot

clean:
	$(RM) *.pdf *.eps *depend.tex *.gpi