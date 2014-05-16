SHELL		= /bin/sh
CC		= gcc
CFLAGS		= -Wall -pedantic -std=c99
LIBRARIES	= -lm -lgsl -lgslcblas
CLEAN_FILES	= kepler.o transits.o *~
DISTCLEAN_FILES	= example *.dat

.PHONY: clean distclean check-syntax

example: example.c kepler.o transits.o
	$(CC) $^ $(CFLAGS) $(LIBRARIES) -o $@

transits.o: transits.c transits.h kepler.o
	$(CC) -c $(CFLAGS) $(LIBRARIES) transits.c

kepler.o: kepler.c kepler.h
	$(CC) -c $(CFLAGS) $(LIBRARIES) kepler.c

clean:
	rm -f $(CLEAN_FILES)

distclean: clean
	rm -f $(DISTCLEAN_FILES)

# per effettuare il controllo della sintassi in Emacs con Flymake
check-syntax:
	$(CC) -o /dev/null $(CFLAGS) -S $(CHK_SOURCES)
