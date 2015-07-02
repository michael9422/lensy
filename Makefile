#
# Note to myself: When using the GNU build tools instead of this simple
# Makefile, I must do these steps:
# 
#	1) Add #include "config.h" to the C source file.
#	2) Create the files "configure.ac", "Makefile.am".
#	3) Run the program "autoreconf --install".
#	4) ./configure
#	5) make
#	6) make distcheck
#
#
all: telescope spectrograph

CC = gcc
INCLUDE = .
CFLAGS = -Wall
PREFIX = /usr/local
CLIBS = -lm -lSDL2

spectrograph: liblensy.a spectrograph.c
	$(CC) spectrograph.c liblensy.a -I$(INCLUDE) $(CFLAGS) $(CLIBS) -o spectrograph

telescope: liblensy.a telescope.c
	$(CC) telescope.c liblensy.a -I$(INCLUDE) $(CFLAGS) $(CLIBS) -o telescope

liblensy.a: lensy.o
	ar crv liblensy.a lensy.o
	ranlib liblensy.a

lensy.o: lensy.c lensy.h list.h
	$(CC) -I$(INCLUDE) $(CFLAGS) -lm -c lensy.c

install: liblensy.a
	@if [ -d $(PREFIX) ]; then \
	   cp liblensy.a $(PREFIX)/lib/; \
	   cp lensy.h $(PREFIX)/include/; \
	   cp list.h $(PREFIX)/include/; \
#	   cp spectrograph $(PREFIX)/sbin/; \
	   cp telescope $(PREFIX)/sbin/; \
	   echo "Installed in $(PREFIX)/sbin"; \
	else \
	   echo "Sorry, $(PREFIX) does not exist"; \
	fi

uninstall:
	-rm $(PREFIX)/lib/liblensy.a
	-rm $(PREFIX)/include/lensy.h
	-rm $(PREFIX)/sbin/telescope
#	-rm $(PREFIX)/sbin/spectrograph

clean:
	-rm liblensy.a
