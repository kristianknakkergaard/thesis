CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: xifilldepend xifilldepend1.tex xifilldepend2.tex

xifilldepend1.tex: xifilldepend Makefile
	echo \
'set terminal epslatex color; 			'\
'set output "$@"; 			   		'\
'if (!exists("MP_LEFT"))   MP_LEFT = .05; 	'\
'if (!exists("MP_RIGHT"))  MP_RIGHT = .95; 	'\
'if (!exists("MP_BOTTOM")) MP_BOTTOM = .1; 	'\
'if (!exists("MP_TOP"))    MP_TOP = .9; 	'\
'if (!exists("MP_GAP"))    MP_GAP = .05; 	'\
'set  multiplot layout 1,2 rowsfirst margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP; '\
'unset grid;  						'\
'set tics out;       			    '\
'unset ylabel;						'\
'set xlabel "$$n$$";				'\
'unset key;							'\
'unset colorbox;					'\
'set xrange  [0:1]; 				'\
'set yrange  [1:7]; 				'\
'set palette defined (0 "white", 1 "red", 2 "blue" );						'\
'set label "$$(i)$$" center at graph 0.2, 0.7142857143 front textcolor rgb "white";			'\
'set label "$$(ii)$$" center at graph 0.5, 0.7142857143 front textcolor rgb "white";			'\
'set ylabel "$$\xi / d$$";													'\
'set title "$$t_2 = t_1, \text{anomalous}$$";								'\
'plot  																		'\
' "$<"  index 4 using 1:2:3 with image pixels;								'\
'set palette defined (0 "white", 1 "red");									'\
'set ytics format " "; 														'\
'unset ylabel; 																'\
'set title "$$t_2 = t_1, \text{normal}$$"; 									'\
'plot  																		'\
' "$<"  index 5 using 1:2:3 with image pixels;								'\
| gnuplot

xifilldepend2.tex: xifilldepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'if (!exists("MP_LEFT"))   MP_LEFT = .05; 	'\
'if (!exists("MP_RIGHT"))  MP_RIGHT = .95; 	'\
'if (!exists("MP_BOTTOM")) MP_BOTTOM = .1; 	'\
'if (!exists("MP_TOP"))    MP_TOP = .9; 	'\
'if (!exists("MP_GAP"))    MP_GAP = .05; 	'\
'set  multiplot layout 1,2 rowsfirst margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP; '\
'unset grid;  						'\
'set tics out;       			    '\
'unset ylabel;						'\
'set xlabel "$$n$$";				'\
'unset key;							'\
'unset colorbox;					'\
'set xrange  [0:1]; 				'\
'set yrange  [1:7]; 				'\
'set palette defined (0 "white", 1 "red", 2 "blue" );						'\
'set ylabel "$$\xi / d$$";													'\
'set title "$$t_2 = t_1, \text{anomalous}$$";								'\
'plot  																		'\
' "$<"  index 0 using 1:2:3 with image pixels;								'\
'set palette defined (0 "white", 1 "red");									'\
'set ytics format " "; 														'\
'unset ylabel; 																'\
'set title "$$t_2 = t_1, \text{normal}$$"; 									'\
'plot  																		'\
' "$<"  index 1 using 1:2:3 with image pixels;								'\
| gnuplot


clean:
	$(RM) *.pdf *.eps *depend.tex *.gpi