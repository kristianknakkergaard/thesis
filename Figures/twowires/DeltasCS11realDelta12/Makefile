CFLAGS = -Wall -Werror -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: ddepend

ddepend: CS11Delta12real
	./$< > ddepend

CS11Delta12real: CS11Delta12real.o

clean:
	$(RM) CS11Delta12real CS11Delta12real.o ddepend