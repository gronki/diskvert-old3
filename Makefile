VERSION  	 	= 170530

prefix 	 	 	= /usr/local
bindir		 	= $(prefix)/bin
datadir	 	 	= $(prefix)/share
includedir 	 	= $(prefix)/include
libdir 	 	 	= $(prefix)/lib
fmoddir			= $(libdir)/finclude

INCLUDE	 	 	= -J. -I. -Ilibconfort
LDFLAGS			= -L. -Llibconfort

FC				= f95
CFLAGS			?= -g -O3 -march=native
FFLAGS			?= $(CFLAGS)

OBJECTS_MATH = bisect.o cgs.o deriv.o eulerintegr.o findzer.o findzer_multi.o histogram.o interpol.o  kramers.o linsect.o random.o rk4integr.o space.o threshold.o
OBJECTS_LIB = $(OBJECTS_MATH) alphadisk.o balance.o globals.o model_m1.o model_ss73.o
OBJECTS_UTIL = results.o summary.o setup.o settings.o

VPATH = src:src/util:src/prog:src/math

PROGRAMS = diskvert-m1 diskvert-ss73 disk-properties
BINARIES = $(addprefix bin/,$(PROGRAMS))

all: $(BINARIES) libdiskvert.so

#################  INSTALACJA  #################

install: all
	install -d $(libdir)
	install libdiskvert.so $(libdir)
	install -d $(libdir)
	install libconfort/libconfort.so $(libdir)
	install -d $(bindir)
	install $(BINARIES) $(bindir)
	install scripts/diskvert-pack $(bindir)

install-user: prefix = $(HOME)/.local
install-user: install
	@echo
	@echo "Program installed in $(bindir) and $(libdir),"
	@echo "Before usage, you must invoke:"
	@echo "   export PATH=\"\$$PATH:$(bindir)\""
	@echo "   export LD_LIBRARY_PATH=\"\$$LD_LIBRARY_PATH:$(libdir)\""
	@echo

#################  PLIKI OBIEKTOW  #################

%.o: %.F90
	$(FC) $(INCLUDE) $(FFLAGS) $(CPPFLAGS) -c -fPIC $< -o $@

$(OBJECTS_LIB) $(OBJECTS_UTIL): kind.o
settings.o: libconfort/libconfort.so
alphadisk.o balance.o: settings.o globals.o
settings.o: globals.o
setup.o: settings.o globals.o
settings.o: results.o
balance.o: globals.o
alphadisk.o: globals.o
random.o: interpol.o
kramers.o: cgs.o
globals.o: $(OBJECTS_MATH)
model_m1.o model_ss73.o: $(OBJECTS_MATH) alphadisk.o balance.o globals.o
diskvert-m1.o: model_m1.o
diskvert-ss73.o: model_ss73.o

#################  PLIKI BINARNE  #################

libdiskvert.so: $(OBJECTS_LIB)
	$(FC) $(LDFLAGS) -shared $^ -o $@

libconfort/libconfort.so:
	$(MAKE) -C libconfort libconfort.so

bin:
	mkdir -p bin

$(BINARIES): bin/%: $(OBJECTS_UTIL) %.o | libconfort/libconfort.so libdiskvert.so bin
	$(FC) $(LDFLAGS) $^ -lconfort -ldiskvert -o $@

#################  SPRZATANIE  #################

clean:
	$(RM) *.mod *.smod *.o diskvert/*.pyc
	$(RM) -r dist *.egg-info build
	$(MAKE) -C libconfort clean
distclean: clean
	$(RM) -r bin *.so
	$(MAKE) -C libconfort distclean
