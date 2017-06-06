VERSION  	 	= 170530

prefix 	 	 	= /usr/local
bindir		 	= $(prefix)/bin
datadir	 	 	= $(prefix)/share
includedir 	 	= $(prefix)/include
libdir 	 	 	= $(prefix)/lib
fmoddir			= $(libdir)/gfortran/modules
pkgconfigdir    = $(libdir)/pkgconfig

INCLUDE	 	 	= -I. -Imod -Ilibconfort
LDFLAGS			= -L. -Llibconfort

FC				= f95
CFLAGS			?= -g -Wall -O3 -march=native -mieee-fp
FFLAGS			?= $(CFLAGS) -Warray-temporaries -std=f2008
override CPPFLAGS += -DVERSION=$(VERSION)

OBJECTS_MATH = bisect.o cgs.o deriv.o eulerintegr.o \
	findzer.o findzer_multi.o histogram.o interpol.o  kramers.o \
	linsect.o random.o rk4integr.o space.o threshold.o
OBJECTS_LIB = alphadisk.o balance.o globals.o model_m1.o \
	model_ss73.o coefficients.o heyney_matrix.o
OBJECTS_UTIL = results.o summary.o setup.o settings.o

VPATH = src:src/util:src/prog:src/math

PROGRAMS = diskvert-m1 diskvert-ss73 disk-properties
BINARIES = $(addprefix bin/,$(PROGRAMS))

all: $(BINARIES) libdiskvert.so

#################  INSTALACJA  #################

install: all
	install -d $(libdir)
	install libdiskvert.so $(libdir)
	install libconfort/libconfort.so $(libdir)
	install -d $(fmoddir)
	install -m 644 mod/*.mod $(fmoddir)
	install -d $(bindir)
	install $(BINARIES) $(bindir)
	install scripts/diskvert-pack $(bindir)
	ldconfig -nv $(libdir)
	install -d $(pkgconfigdir)
	@echo "Name: diskvert" | tee $(pkgconfigdir)/diskvert.pc
	@echo "Version: $(VERSION)" | tee -a $(pkgconfigdir)/diskvert.pc
	@echo "Description: compute vertical structure of accretion disks" | tee -a $(pkgconfigdir)/diskvert.pc
	@echo "Cflags: -I$(fmoddir)" | tee -a $(pkgconfigdir)/diskvert.pc
	@echo "Libs: -L$(libdir) -ldiskvert" | tee -a $(pkgconfigdir)/diskvert.pc
	@echo "export PATH=\"\$$PATH:$(bindir)\"" | tee diskvert-paths.sh
	@echo "export LD_LIBRARY_PATH=\"\$$LD_LIBRARY_PATH:$(libdir)\"" | tee -a diskvert-paths.sh
	@echo "export PKG_CONFIG_PATH=\"\$$PKG_CONFIG_PATH:$(pkgconfigdir)\"" | tee -a diskvert-paths.sh
	@echo
	@echo "Program installed in $(bindir) and $(libdir). Invoke one of the following:"
	@echo
	@echo "# activate for the current session"
	@echo ". diskvert-paths.sh"
	@echo "# activate for the current user"
	@echo "cat diskvert-paths.sh | tee -a ~/.bashrc"
	@echo "# activate globally for the whole system"
	@echo "cat diskvert-paths.sh | sudo tee -a /etc/profile.d/diskvert-paths.sh"
	@echo
	@echo "You might want to install Python scripts as well:"
	@echo "python setup.py install --user"
	@echo

install-user: prefix = $(HOME)/.local
install-user: install


#################  PLIKI OBIEKTOW  #################

%.o: %.F90
	$(FC) $(INCLUDE) $(FFLAGS) $(CPPFLAGS) -c -fPIC $< -o $@
$(OBJECTS_LIB) precision.o: %.o: %.F90 | mod
	$(FC) $(INCLUDE) -Jmod $(FFLAGS) $(CPPFLAGS) -c -fPIC $< -o $@

$(OBJECTS_MATH) $(OBJECTS_LIB) $(OBJECTS_UTIL): precision.o
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
coefficients.o: globals.o cgs.o src/coefficients.F90
heyney_matrix.o: globals.o cgs.o coefficients.o

src/coefficients.F90: $(wildcard generate_coefficients/*.py)
	python generate_coefficients

#################  PLIKI BINARNE  #################

libdiskvert.so: $(OBJECTS_MATH) $(OBJECTS_LIB) precision.o
	$(FC) $(LDFLAGS) -shared $^ -o $@

libconfort/libconfort.so:
	$(MAKE) -C libconfort libconfort.so

bin mod:
	mkdir -p $@

$(BINARIES): bin/%: $(OBJECTS_UTIL) %.o | libconfort/libconfort.so libdiskvert.so bin
	$(FC) $(LDFLAGS) $^ -lconfort -ldiskvert -o $@

#################  SPRZATANIE  #################

clean:
	$(RM) *.mod *.smod *.o diskvert/*.pyc generate_coefficients/*.pyc
	$(RM) -r dist *.egg-info build
	$(MAKE) -C libconfort clean
distclean: clean
	$(RM) -r bin mod
	$(RM) *.so diskvert-paths.sh
	$(MAKE) -C libconfort distclean
