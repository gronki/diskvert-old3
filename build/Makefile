VERSION = 180613

#------------------------------------------------------------------------------#
#----------------------------- INSTALLATION PATHS -----------------------------#
#------------------------------------------------------------------------------#

prefix = /usr/local
bindir = $(prefix)/bin
libdir = $(prefix)/lib
includedir = $(prefix)/include
fmoddir = $(includedir)

#------------------------------------------------------------------------------#

override INCLUDE += -I../libconfort
override LDFLAGS += -L.

#------------------------------------------------------------------------------#

BUILD = release
VENDOR = gnu

#------------------------------------------------------------------------------#
#---------------------------- CFLAGS AND COMPILERS ----------------------------#
#------------------------------------------------------------------------------#
ifeq ($(VENDOR), intel)
#------------------------------------------------------------------------------#
# Intel compilers

CC := icc
FC := ifort

ifeq ($(BUILD), debug)
FFLAGS := -g -Og -C -fp-model except -warn all
else
FFLAGS := -O3 -xhost -ipo -warn
endif

#------------------------------------------------------------------------------#
else ifeq ($(VENDOR), pgi)
#------------------------------------------------------------------------------#
# PGI compilers

$(error PGI compilers are not supported for now)

#------------------------------------------------------------------------------#
else
#------------------------------------------------------------------------------#
# GNU compilers (default)

CC := gcc
FC := gfortran

ifeq ($(BUILD), debug)
FFLAGS := -g -Og -fcheck=all -fbacktrace
else
FFLAGS ?= -march=native
FFLAGS += -O3 -flto
endif

FFLAGS += -fimplicit-none -ffree-line-length-none -pedantic -fno-backtrace
FFLAGS += -Wall -Wno-unused-dummy-argument -Warray-temporaries -Wrealloc-lhs

#------------------------------------------------------------------------------#
endif
#------------------------------------------------------------------------------#
#------------------------- CFLAGS AND COMPILERS - END -------------------------#
#------------------------------------------------------------------------------#

export CC FC

#------------------------------------------------------------------------------#
#-------------------------- SOURCES AND OBJECT FILES --------------------------#
#------------------------------------------------------------------------------#

# search for sources in these directories
VPATH := ../src:../src/util:../src/prog:../src/math

# objects that go into the shared library
objects := balance.o globals.o grid.o heatbil.o ss73solution.o modelmag.o \
	alphadisk.o relaxation.o rk4settings.o
# various math utilities
objects += cgs.o threshold.o rk4integr.o eulerintegr.o ramps.o integrate.o \
	interpol.o deriv.o

objects_lib := bisect.o findzer.o findzer_multi.o histogram.o linsect.o \
	random.o relaxation_c.o
objects_util := results.o settings.o fileunits.o

#------------------------------------------------------------------------------#

# if using intel compiler, use the MKL instead of standard LAPACK
ifeq ($(VENDOR), intel)
LDLIBS += -mkl=sequential
else
LDLIBS += -llapack
endif

#------------------------------------------------------------------------------#

programs = dv-alpha dv-alpha-rx dv-mag dv-mag-rx dv-rad1

#------------------------------------------------------------------------------#
#-------------------------------- BUILD RULES ---------------------------------#
#------------------------------------------------------------------------------#

all: $(programs)

lib: override FFLAGS += -fPIC
lib: libdiskvert.so

#------------------------------------------------------------------------------#

install: all
	install -d $(DESTDIR)$(bindir)
	install $(programs) $(DESTDIR)$(bindir)
	ln -sf $(bindir)/dv-mag-rx $(DESTDIR)$(bindir)/diskvert

install-lib: lib
	install -d $(DESTDIR)$(libdir)
	install libdiskvert.so $(DESTDIR)$(libdir)
	install -d $(DESTDIR)$(fmoddir)/diskvert
	install -m 644 slf_{cgs,rk4integr,threshold}.mod \
	 	{globals,settings,grid,rk4settings}.mod \
		{alphadisk,modelmag,relaxation,ss73solution}.mod \
		$(DESTDIR)$(fmoddir)/diskvert

#------------------------------------------------------------------------------#

%.o: %.F90
	$(FC) $(INCLUDE) $(CPPFLAGS) $(FFLAGS) -c $< -o $@
%.o: %.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c $< -o $@

include deps.inc

relaxation.o    : mrxcoeff.fi mrxdims.fi mrxhash.fi mrxptrs.fi
settings.o      : ../libconfort/libconfort.a

#------------------------------------------------------------------------------#

$(programs): %: $(objects) $(objects_util) %.o ../libconfort/libconfort.a
	$(FC) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

#------------------------------------------------------------------------------#

libdiskvert.so: $(objects) $(objects_lib)
	$(FC) $(FFLAGS) $(LDFLAGS) -shared $^ $(LDLIBS) -o $@

#------------------------------------------------------------------------------#

../libconfort/libconfort.a:
	$(MAKE) -C ../libconfort CC=$(CC) FC=$(FC) libconfort.a

#------------------------------------------------------------------------------#

clean:
	$(RM) *.mod *.smod *.a *.o
	$(RM) $(programs) *.so
	$(MAKE) -C ../libconfort distclean
