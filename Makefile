VERSION  	 	= 170623

prefix 	 	 	= /usr/local
bindir		 	= $(prefix)/bin
datadir	 	 	= $(prefix)/share
includedir 	 	= $(prefix)/include
libdir 	 	 	= $(prefix)/lib
fmoddir			= $(libdir)/gfortran/modules
pkgconfigdir    = $(libdir)/pkgconfig

INCLUDE	 	 	= -Ilibconfort
LDFLAGS			= -L.
LDLIBS			= -lopenblas

FC				= f95
CFLAGS			?= -g -Wall -O3 -march=native -mieee-fp
FFLAGS			?= $(CFLAGS) -Warray-temporaries -Wpedantic
override FFLAGS += -fexternal-blas
override CPPFLAGS += -DVERSION=$(VERSION)

OBJECTS_LAPACK = $(addsuffix .o,$(basename $(notdir \
	$(wildcard src/lapack/*.[fF]))))

OBJECTS_MATH = $(addsuffix .o,$(basename $(notdir \
	$(wildcard src/math/*.[fF]90))))
OBJECTS_LIB =  $(addsuffix .o,$(basename $(notdir \
	$(wildcard src/*.[fF]90))))
OBJECTS_UTIL = $(addsuffix .o,$(basename $(notdir \
	$(wildcard src/util/*.[fF]90))))

VPATH = src:src/util:src/prog:src/math:src/lapack

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

include make_dependencies.inc

relaxation.o    : lapack.a
settings.o      : libconfort.a
src/coefficients.f90: $(wildcard generate_coefficients/*.py)
	python generate_coefficients

#################  PLIKI BINARNE  #################

libdiskvert.so: $(OBJECTS_MATH) $(OBJECTS_LIB) lapack.a
	$(FC) $(LDFLAGS) -shared $^ $(LDLIBS) -o $@

lapack.a: $(OBJECTS_LAPACK)
	$(AR) rcs $@ $^
.INTERMEDIATE: $(OBJECTS_LAPACK)

libconfort.a:
	$(MAKE) -C libconfort libconfort.a
	ln -s libconfort/libconfort.a $@

bin:
	mkdir -p $@

$(BINARIES): bin/%: $(OBJECTS_UTIL) %.o libconfort.a | libdiskvert.so bin
	$(FC) $(LDFLAGS) $^ -ldiskvert $(LDLIBS) -o $@

#################  SPRZATANIE  #################

clean:
	$(RM) *.mod *.smod *.a *.o diskvert/*.pyc generate_coefficients/*.pyc
	$(RM) -r dist *.egg-info build
	$(MAKE) -C libconfort clean
distclean: clean
	$(RM) -r bin mod
	$(RM) *.so diskvert-paths.sh
	$(MAKE) -C libconfort distclean
