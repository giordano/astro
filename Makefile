SHELL		= /bin/sh
CC		= gcc-4.4
CFLAGS		= -Wall -pedantic
LIBRARIES	= -lm -lgsl -lgslcblas
CLEAN_FILES	= keplero.o transiti.o *~
DISTCLEAN_FILES	= esempio *.dat

.PHONY: clean distclean check-syntax

esempio: esempio.c keplero.o transiti.o
	$(CC) $(CFLAGS) $(LIBRARIES) -o $@ $^

transiti.o: transiti.c transiti.h keplero.o
	$(CC) -c $(CFLAGS) $(LIBRARIES) transiti.c

keplero.o: keplero.c keplero.h
	$(CC) -c $(CFLAGS) $(LIBRARIES) keplero.c

clean:
	rm -f $(CLEAN_FILES)

distclean: clean
	rm -f $(DISTCLEAN_FILES)

# per effettuare il controllo della sintassi in Emacs con Flymake
check-syntax:
	$(CC) -o /dev/null $(CFLAGS) -S $(CHK_SOURCES)
