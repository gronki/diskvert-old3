VERSION  	 	= 170623

prefix 	 	 	= /usr/local
bindir		 	= $(prefix)/bin
datadir	 	 	= $(prefix)/share
includedir 	= $(prefix)/include
libdir 	 	 	= $(prefix)/lib
fmoddir			= $(libdir)/gfortran/modules
pkgconfigdir= $(libdir)/pkgconfig

INCLUDE	 	 	= -Ilibconfort
LDFLAGS			= -L.
LDLIBS			= -lopenblas

FC					= f95
CFLAGS			?= -g -Wall -O3 -mieee-fp
FFLAGS			?= $(CFLAGS) -Warray-temporaries -Wpedantic \
			-Wno-unused-dummy-argument
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

coeffincludes = coefficients.fi mrxdims.fi mrxhash.fi mrxname.fi mrxptrs.fi

VPATH = src:src/util:src/prog:src/math:src/lapack

PROGRAMS = $(basename $(notdir $(wildcard src/prog/*.[fF]90)))
BINARIES = $(addprefix bin/,$(PROGRAMS))

all: $(BINARIES) libdiskvert.so

#################  INSTALACJA  #################

install: all
	install -d $(libdir)
	install libdiskvert.so $(libdir)
	install -d $(fmoddir)/diskvert
	install -m 644 slf_{cgs,rk4integr,threshold}.mod \
	 	{globals,settings,fileunits,grid,rk4settings}.mod \
		{alphadisk,alphasimp,modelmag,relaxation,ss73solution}.mod \
		$(fmoddir)/diskvert
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

$(coeffincludes): relaxation-coefficients.py
	python relaxation-coefficients.py

relaxation.o: $(coeffincludes)

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
