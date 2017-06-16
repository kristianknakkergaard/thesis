CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: kdepend.tex kdepend1.tex  kdepend2.tex kdepend dispersion

kdepend.tex: kdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'if (!exists("MP_LEFT"))   MP_LEFT = .05; 	'\
'if (!exists("MP_RIGHT"))  MP_RIGHT = .95; 	'\
'if (!exists("MP_BOTTOM")) MP_BOTTOM = .15; '\
'if (!exists("MP_TOP"))    MP_TOP = .9; 	'\
'if (!exists("MP_GAP"))    MP_GAP = .1; 	'\
'set  multiplot layout 1,2 rowsfirst margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP; '\
'set grid lw 2;  															'\
'set tics out;       			    										'\
'set key top right Left;													'\
'set key opaque width -10;													'\
'set key spacing 1.5; 														'\
'set arrow from 0,0 to 0,1.5 nohead lc rgb "black"; 							'\
'set arrow from -1.5,0 to 1.5,0 nohead lc rgb "black"; 							'\
'unset ylabel;																'\
'set xlabel "$$k / k_F$$";													'\
'set xtics -1.5, 0.5, 1.5; 														'\
'set ylabel "$$E / \epsilon_{F,0}$$";													'\
'set title "$$Imaginary interwire pairing $$"; 											'\
'plot [-1.5:1.5][0:1.5]																		'\
' "$<" index 0    using 1:2 title "$$|\varepsilon_k|$$" with lines lw 2 lc rgb "black"	'\
',"$<" index 0    using 1:3 title "$$E_{F,k}$$" with lines lw 2 lc rgb "red";    			'\
'set ytics format " "; 														'\
'unset ylabel; 																'\
'set title "$$Real interwire pairing $$"; 												'\
'plot [-1.5:1.5][0:1.5]																		'\
' "$<" index 0    using 1:2 title "$$|\varepsilon_k|$$" with lines lw 2 lc rgb "black"	'\
',"$<" index 0    using 1:4 title "$$E^{+}_{F,k}$$" with lines lw 2 lc rgb "red"				'\
',"$<" index 0    using 1:5 title "$$E^{-}_{F,k}$$" with lines lw 2 lc rgb "blue";			'\
| gnuplot

kdepend1.tex: kdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid lw 2;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$E_{F,k}/epsilon_{F,0}$$";'\
'set key top right Left;					'\
'set key opaque width -1;						'\
'set arrow from 0,0 to 0,3 nohead lc rgb "black"; 	'\
'set arrow from -2,0 to 2,0 nohead lc rgb "black"; 	'\
'set title "$$Imaginary interwire pairing $$"; 			'\
'plot [-2:2][0:3]			'\
' "$<" index 0    using 1:2 title "$$|\varepsilon_k|$$" with lines lw 2 lc rgb "black"'\
',"$<" index 0    using 1:3 title "$$E_{F,k}$$" with lines lw 2 lc rgb "red"'\
| gnuplot

kdepend2.tex: kdepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'set grid;       			   		'\
'set tics out;       			    '\
'set xlabel "$$k/k_F$$"; 		    '\
'set ylabel "$$E_{F,k}/epsilon_{F,0}$$";'\
'set key top right Left;					'\
'set key opaque width -1;						'\
'set arrow from 0,0 to 0,3 nohead lc rgb "black"; 	'\
'set arrow from -2,0 to 2,0 nohead lc rgb "black"; 	'\
'set title "$$Real interwire pairing $$"; 			'\
'plot [-2:2][0:3]			'\
' "$<" index 0    using 1:2 title "$$|\varepsilon_k|$$" with lines lw 2 lc rgb "black"'\
',"$<" index 0    using 1:4 title "$$E^{+}_{F,k}$$" with lines lw 2 lc rgb "red"'\
',"$<" index 0    using 1:5 title "$$E^{-}_{F,k}$$" with lines lw 2 lc rgb "blue"'\
| gnuplot


kdepend: dispersion
	./$< 1> kdepend


dispersion: dispersion.o

clean:
	$(RM) dispersion dispersion.o *.pdf *.eps *depend.tex kdepend *.gpi