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
'if (!exists("MP_BOTTOM")) MP_BOTTOM = .1; 	'\
'if (!exists("MP_TOP"))    MP_TOP = .9; 	'\
'if (!exists("MP_GAP"))    MP_GAP = .05; 	'\
'set  multiplot layout 1,2 rowsfirst margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP; '\
'unset grid;       			   												'\
'set tics out;       			   											'\
'set palette defined (0 "white", 1 "red" );									'\
'set ylabel "$$\xi / d$$";													'\
'set xlabel "$$n$$";														'\
'unset key;																	'\
'unset colorbox;															'\
'set title "$$G = 4$$";														'\
'plot  [0:1][1:7]															'\
' "$<"  index 0 using 1:2:3 with image pixels;								'\
'unset ylabel;																'\
'set ytics format " "; 														'\
'set title "$$nu = 8$$"; 													'\
'plot  [0:1][1:7]															'\
' "$<"  index 1 using 1:2:3 with image pixels;								'\
| gnuplot

clean:
	$(RM) *.pdf *.eps *depend.tex *.gpi