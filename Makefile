SHELL		= /bin/sh
HEADERS		= keplero.h transiti.h
CFLAGS		= -Wall -pedantic
LIBRARIES	= -lm -lgsl -lgslcblas
CLEAN_FILES	= keplero.o transiti.o *~
DISTCLEAN_FILES	= $(EXES) *.dat

.PHONY: clean distclean check-syntax

esempio: esempio.c keplero.o transiti.o
	gcc $(CFLAGS) $(LIBRARIES) -o $@ $^

transiti.o: transiti.c transiti.h keplero.o
	gcc -c $(CFLAGS) $(LIBRARIES) transiti.c

keplero.o: keplero.c keplero.h
	gcc -c $(CFLAGS) $(LIBRARIES) keplero.c

clean:
	rm -f $(CLEAN_FILES)

distclean: clean
	rm -f $(DISTCLEAN_FILES)

# per effettuare il controllo della sintassi in Emacs con Flymake
check-syntax:
	gcc -o /dev/null $(CFLAGS) -S $(CHK_SOURCES)
