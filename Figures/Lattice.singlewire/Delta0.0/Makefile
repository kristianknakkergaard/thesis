CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: xifilldepend.tex

xifilldepend.tex: xifilldepend Makefile
	echo \
'set terminal epslatex; 			'\
'set output "$@"; 			   		'\
'if (!exists("MP_LEFT"))   MP_LEFT = .05; 	'\
'if (!exists("MP_RIGHT"))  MP_RIGHT = .95; 	'\
'if (!exists("MP_BOTTOM")) MP_BOTTOM = .15; '\
'if (!exists("MP_TOP"))    MP_TOP = .9; 	'\
'if (!exists("MP_GAP"))    MP_GAP = .05; 	'\
'set  multiplot layout 1,2 rowsfirst margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP; '\
'unset grid;       			   												'\
'set tics out;       			   											'\
'set palette defined (0 "white", 1 "red" );									'\
'set label 1 "$$nu = 0$$" rotate 	center at graph 0.07, 0.7142857143 front; 	'\
'set label 2 "$$nu = 0$$" rotate 	center at graph 0.93, 0.7142857143 front; 	'\
'set label 3 "$$nu = 1$$" 	 		center at graph 0.50, 0.7142857143 front; 	'\
'set ylabel "$$\xi / d$$";													'\
'set xlabel "$$n$$";														'\
'unset key;																	'\
'unset colorbox;															'\
'set title "$$G = 4$$";														'\
'plot  [0:1][1:7]															'\
' "$<"  index 0 using 1:2:3 with image pixels;								'\
'unset ylabel;																'\
'set label 1 "$$nu = 0$$" 		 	center at graph 0.105, 0.7142857143 front; 	'\
'set label 2 "$$nu = 0$$"  			center at graph 0.895, 0.7142857143 front; 	'\
'set label 3 "$$nu = 1$$" 	 		center at graph 0.500, 0.7142857143 front; 	'\
'set ytics format " "; 														'\
'set title "$$G = 8$$"; 													'\
'plot  [0:1][1:7]															'\
' "$<"  index 1 using 1:2:3 with image pixels;								'\
| gnuplot

clean:
	$(RM) *.pdf *.eps *depend.tex *.gpi