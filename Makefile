SHELL		= /bin/sh
CC		= gcc
CFLAGS		= -Wall -pedantic -std=c99
AR		= ar
LIBRARIES	= -lm -lgsl -lgslcblas
CLEAN_FILES	= *~
DISTCLEAN_FILES	= example *.dat *.a *.o

.PHONY: all install uninstall clean distclean check-syntax

all: libgastro.a libgastro-fortran.a

libgastro.a: kepler.o transits.o
	$(AR) rcs $@ $^

libgastro-fortran.a: kepler-fortran.o transits-fortran.o
	$(AR) rcs $@ $^

example: example.c kepler.o transits.o
	$(CC) $^ $(CFLAGS) $(LIBRARIES) -o $@

transits.o: transits.c transits.h kepler.o
	$(CC) -c $(CFLAGS) $(LIBRARIES) transits.c

transits-fortran.o: transits-fortran.c transits-fortran.h transits.o
	$(CC) -c $(CFLAGS) $(LIBRARIES) transits-fortran.c

kepler.o: kepler.c kepler.h
	$(CC) -c $(CFLAGS) $(LIBRARIES) kepler.c

kepler-fortran.o: kepler-fortran.c kepler-fortran.h kepler.o
	$(CC) -c $(CFLAGS) $(LIBRARIES) kepler-fortran.c

install: libgastro.a libgastro-fortran.a kepler.h kepler-fortran.h \
		transits.h transits-fortran.h lensing.f
	install -m 644 libgastro.a /usr/local/lib
	install -m 644 libgastro-fortran.a /usr/local/lib
	install -m 644 kepler.h /usr/local/include
	install -m 644 kepler-fortran.h /usr/local/include
	install -m 644 transits.h /usr/local/include
	install -m 644 transits-fortran.h /usr/local/include
	install -m 644 lensing.f /usr/local/include

uninstall:
	rm -rf /usr/local/lib/libgastro.a /usr/local/lib/libgastro-fortran.a \
	/usr/local/include/kepler.h /usr/local/include/kepler-fortran.h \
	/usr/local/include/transits.h /usr/local/include/transits-fortran.h

clean:
	rm -f $(CLEAN_FILES)

distclean: clean
	rm -f $(DISTCLEAN_FILES)

check-syntax:
	$(CC) -o /dev/null $(CFLAGS) -S $(CHK_SOURCES)
