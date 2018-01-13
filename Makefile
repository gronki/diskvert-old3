VERSION = 180113

#------------------------------------------------------------------------------#

prefix = /usr/local
bindir = $(prefix)/bin
datadir = $(prefix)/share
includedir = $(prefix)/include
libdir = $(prefix)/lib
fmoddir = $(libdir)/gfortran/modules
pkgconfigdir = $(libdir)/pkgconfig

#------------------------------------------------------------------------------#

INCLUDE = -Ilibconfort
LDFLAGS = -L.

#------------------------------------------------------------------------------#

FC = gfortran
FFLAGS ?= -g -O2

ifeq ($(FC),gfortran)
FFLAGS += -fno-unsafe-math-optimizations -Wall -Warray-temporaries -Wrealloc-lhs-all
endif

ifeq ($(FC),ifort)
FFLAGS += -fp-model precise -Winline
endif

#------------------------------------------------------------------------------#

# search for sources in these directories
VPATH = src:src/util:src/prog:src/math

# objects that go into the shared library
OBJECTS =  $(addsuffix .o,$(basename $(notdir \
			$(wildcard src/*.[fF]90))))
# various math utilities
OBJECTS += $(addsuffix .o,$(basename $(notdir \
			$(wildcard src/math/*.[fF]90))))
# objects that are needed for binary programs only
OBJECTS_UTIL = $(addsuffix .o,$(basename $(notdir \
			$(wildcard src/util/*.[fF]90))))

#------------------------------------------------------------------------------#

# if using intel compiler, use the MKL instead of standard LAPACK
ifeq ($(FC),ifort)
LDLIBS += -mkl=sequential
else
# openblas is much faster than standard BLAS and allows
# parallel processing, just like ifort and MKL
LDLIBS += -lopenblas
# compile lapack from source. for unknown reasons, this works faster on PSK
OBJECTS += $(addsuffix .o,$(basename $(notdir \
			$(wildcard src/lapack/*.[fF]))))
VPATH := $(VPATH):src/lapack
endif

#------------------------------------------------------------------------------#

PROGRAMS = $(basename $(notdir $(wildcard src/prog/*.[fF]90)))
BINARIES = $(addprefix bin/,$(PROGRAMS))

#------------------------------------------------------------------------------#

all: $(BINARIES) libdiskvert.so

#------------------------------------------------------------------------------#

install: all
	install -d $(libdir)
	install libdiskvert.so $(libdir)
	install -d $(fmoddir)/diskvert
	install -m 644 slf_{cgs,rk4integr,threshold}.mod \
	 	{globals,settings,fileunits,grid,rk4settings}.mod \
		{alphadisk,modelmag,relaxation,ss73solution}.mod \
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
	@echo "export PATH=\"$(bindir):\$$PATH\"" | tee diskvert-paths.sh
	@echo "export LD_LIBRARY_PATH=\"$(libdir):\$$LD_LIBRARY_PATH\"" | tee -a diskvert-paths.sh
	@echo "export PKG_CONFIG_PATH=\"$(pkgconfigdir):\$$PKG_CONFIG_PATH\"" | tee -a diskvert-paths.sh
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

#------------------------------------------------------------------------------#

%.o %.mod: %.F90
	$(FC) $(INCLUDE) $(FFLAGS) $(CPPFLAGS) -c -fPIC $< -o $@
%.o %.mod: %.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c -fPIC $< -o $@
%.o: %.F
	$(FC) $(INCLUDE) $(FFLAGS) $(CPPFLAGS) -c -fPIC $< -o $@
%.o: %.f
	$(FC) $(INCLUDE) $(FFLAGS) -c -fPIC $< -o $@

include make_dependencies.inc

relaxation.o    : mrxcoeff.fi mrxdims.fi mrxhash.fi mrxptrs.fi
settings.o      : libconfort.a

#------------------------------------------------------------------------------#

libdiskvert.so: $(OBJECTS)
	$(FC) $(LDFLAGS) -shared $^ $(LDLIBS) -o $@

libconfort.a:
	$(MAKE) -C libconfort libconfort.a
	ln -sf libconfort/libconfort.a $@

#------------------------------------------------------------------------------#

bin:
	mkdir -p $@

$(BINARIES): bin/%: $(OBJECTS_UTIL) %.o libconfort.a | libdiskvert.so bin
	$(FC) $(LDFLAGS) $^ -ldiskvert $(LDLIBS) -o $@

#------------------------------------------------------------------------------#

dist: distclean
	tar czfv pydiskvert-$(VERSION).tar.gz 				\
		$(shell git ls-files --exclude-standard)  	\
		--transform "s/^/pydiskvert-$(VERSION)\//"

clean:
	$(RM) *.mod *.smod *.a *.o diskvert/*.pyc generate_coefficients/*.pyc
	$(MAKE) -C libconfort clean
distclean: clean
	$(RM) -r dist *.egg-info build
	$(RM) -r bin mod
	$(RM) *.so diskvert-paths.sh
	$(MAKE) -C libconfort distclean
