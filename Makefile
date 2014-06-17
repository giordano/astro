SHELL		= /bin/sh
CC		= gcc
CFLAGS		= -Wall -Wextra -pedantic -std=c99
AR		= ar
LIBRARIES	= -lm -lgsl -lgslcblas
CLEAN_FILES	= *~
DISTCLEAN_FILES	= example *.dat *.a *.o
INSTALL_LIBS	= libgastro.a
LIBS_DIR	= /usr/local/lib
INSTALLED_LIBS	= $(patsubst %, $(LIBS_DIR)/%, $(INSTALL_LIBS))
INSTALL_HEADERS	= kepler.h kepler-fortran.h transits.h transits-fortran.h \
		  lensing.h lensing-fortran.h lensing.f
HEADERS_DIR	= /usr/local/include
INSTALLED_HEADERS = $(patsubst %, $(HEADERS_DIR)/%, $(INSTALL_HEADERS))
INSTALLED_FILES = $(INSTALLED_LIBS) $(INSTALLED_HEADERS)

.PHONY: all install uninstall clean distclean check-syntax

all: libgastro.a

libgastro.a: kepler.o transits.o lensing.o kepler-fortran.o transits-fortran.o \
	     lensing-fortran.o
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

lensing.o: lensing.c lensing.h
	$(CC) -c $(CFLAGS) $(LIBRARIES) lensing.c

lensing-fortran.o: lensing-fortran.c lensing-fortran.h lensing.o
	$(CC) -c $(CFLAGS) $(LIBRARIES) lensing-fortran.c

$(LIBS_DIR)/%: %
	install -m 644 $^ $(LIBS_DIR)

$(HEADERS_DIR)/%: %
	install -m 644 $^ $(HEADERS_DIR)

install: $(INSTALLED_FILES)

uninstall:
	rm -rf $(INSTALLED_FILES)

clean:
	rm -f $(CLEAN_FILES)

distclean: clean
	rm -f $(DISTCLEAN_FILES)

check-syntax:
	$(CC) -o /dev/null $(CFLAGS) -S $(CHK_SOURCES)
