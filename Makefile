VERSION  	 	= 170530

prefix 	 	 	= /usr/local
bindir		 	= $(prefix)/bin
datadir	 	 	= $(prefix)/share
includedir 	 	= $(prefix)/include
libdir 	 	 	= $(prefix)/lib
fmoddir			= $(libdir)/gfortran/modules
pkgconfigdir    = $(libdir)/pkgconfig

INCLUDE	 	 	= -I. -Ilibconfort
LDFLAGS			= -L. -Llibconfort

FC				= f95
CFLAGS			?= -g -Wall -O3 -march=native -mieee-fp
FFLAGS			?= $(CFLAGS) -Warray-temporaries -fexternal-blas
override CPPFLAGS += -DVERSION=$(VERSION)

OBJECTS_BLAS = dasum.o daxpy.o dcabs1.o dcopy.o ddot.o dgbmv.o dgemm.o dgemv.o dger.o \
	dnrm2.o drot.o drotg.o drotm.o drotmg.o dsbmv.o dscal.o dsdot.o dspmv.o \
	dspr2.o dspr.o dswap.o dsymm.o dsymv.o dsyr2.o dsyr2k.o dsyr.o dsyrk.o \
	dtbmv.o dtbsv.o dtpmv.o dtpsv.o dtrmm.o dtrmv.o dtrsm.o dtrsv.o \
	dzasum.o dznrm2.o icamax.o idamax.o isamax.o izamax.o lsame.o sasum.o \
	saxpy.o scabs1.o scasum.o scnrm2.o scopy.o sdot.o sdsdot.o sgbmv.o \
	sgemm.o sgemv.o sger.o snrm2.o srot.o srotg.o srotm.o srotmg.o ssbmv.o \
	sscal.o sspmv.o sspr2.o sspr.o sswap.o ssymm.o ssymv.o ssyr2.o ssyr2k.o \
	ssyr.o ssyrk.o stbmv.o stbsv.o stpmv.o stpsv.o strmm.o strmv.o strsm.o \
	strsv.o xerbla_array.o xerbla.o

OBJECTS_LAPACK = dgbsv.o dgbtf2.o dgbtrf.o dgbtrs.o dgesv.o dgetf2.o dgetrf2.o \
	dgetrf.o dgetrs.o dlamch.o dlaswp.o ieeeck.o ilaenv.o iparam2stage.o iparmq.o \
	lsame.o xerbla.o $(OBJECTS_BLAS)

OBJECTS_MATH = bisect.o cgs.o deriv.o eulerintegr.o \
	findzer.o findzer_multi.o histogram.o interpol.o  kramers.o \
	linsect.o random.o rk4integr.o space.o threshold.o

OBJECTS_LIB = alphadisk.o balance.o globals.o model_m1.o \
	model_ss73.o coefficients.o relaxation.o

OBJECTS_UTIL = results.o summary.o setup.o settings.o

VPATH = src:src/util:src/prog:src/math:src/lapack:src/lapack/blas

PROGRAMS = diskvert-m1 diskvert-ss73 disk-properties
BINARIES = $(addprefix bin/,$(PROGRAMS))

all: $(BINARIES) libdiskvert.so

#################  INSTALACJA  #################

install: all
	install -d $(libdir)
	install libdiskvert.so $(libdir)
	install -d $(fmoddir)/diskvert
	install -m 644 slf_{kramers,cgs,bisect,space,threshold,rk4integr}.mod \
	 	{alphadisk,model_m1,model_ss73,relaxation,relax_coefficients}.mod \
		{globals,precision,settings,setup}.mod $(fmoddir)/diskvert
	install -d $(bindir)
	install $(BINARIES) $(bindir)
	install scripts/diskvert-pack $(bindir)
	ldconfig -nv $(libdir)
	install -d $(pkgconfigdir)
	@echo "Name: diskvert" | tee $(pkgconfigdir)/diskvert.pc
	@echo "Version: $(VERSION)" | tee -a $(pkgconfigdir)/diskvert.pc
	@echo "Description: compute vertical structure of accretion disks" | tee -a $(pkgconfigdir)/diskvert.pc
	@echo "Cflags: -I$(fmoddir)/diskvert" | tee -a $(pkgconfigdir)/diskvert.pc
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
%.o: %.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c -fPIC $< -o $@
%.o: %.F
	$(FC) $(INCLUDE) $(FFLAGS) $(CPPFLAGS) -c -fPIC $< -o $@
%.o: %.f
	$(FC) $(INCLUDE) $(FFLAGS) -c -fPIC $< -o $@

$(OBJECTS_MATH) $(OBJECTS_LIB) $(OBJECTS_UTIL): precision.o
settings.o: libconfort.a
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
relaxation.o: globals.o cgs.o coefficients.o $(OBJECTS_LAPACK)

src/coefficients.F90: $(wildcard generate_coefficients/*.py)
	python generate_coefficients

#################  PLIKI BINARNE  #################

libdiskvert.so: $(OBJECTS_MATH) $(OBJECTS_LIB) $(OBJECTS_LAPACK) precision.o
	$(FC) $(LDFLAGS) -shared $^ -o $@

libconfort.a:
	$(MAKE) -C libconfort libconfort.a
	ln -s libconfort/libconfort.a $@

bin:
	mkdir -p $@

$(BINARIES): bin/%: $(OBJECTS_UTIL) %.o libconfort.a | libdiskvert.so bin
	$(FC) $(LDFLAGS) $^ -ldiskvert -o $@

#################  SPRZATANIE  #################

clean:
	$(RM) *.mod *.smod *.a *.o diskvert/*.pyc generate_coefficients/*.pyc
	$(RM) -r dist *.egg-info build
	$(MAKE) -C libconfort clean
distclean: clean
	$(RM) -r bin mod
	$(RM) *.so diskvert-paths.sh
	$(MAKE) -C libconfort distclean
