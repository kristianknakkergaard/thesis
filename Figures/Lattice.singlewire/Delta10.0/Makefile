CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: xifilldepend Delta

xifilldepend: Delta
	./$< 1> xifilldepend

Delta: Delta.o

clean:
	$(RM) Delta Delta.o *.pdf *.eps *depend.tex xifilldepend *.gpi







 