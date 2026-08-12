CC ?= gcc
CFLAGS ?= -Wall
LDLIBS = -lm

PROGRAMS = locus matrix vis-spec wl

.PHONY: all clean

all: $(PROGRAMS)

locus: src/locus.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

matrix: src/matrix.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

vis-spec: src/vis-spec.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

wl: src/wl.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

clean:
	rm -f $(PROGRAMS) out.bmp
