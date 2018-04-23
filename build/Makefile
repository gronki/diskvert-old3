VERSION = 180302

#------------------------------------------------------------------------------#

prefix = /usr/local
bindir = $(prefix)/bin
datadir = $(prefix)/share
includedir = $(prefix)/include
libdir = $(prefix)/lib
fmoddir = $(libdir)/gfortran/modules

#------------------------------------------------------------------------------#

override INCLUDE += -I../libconfort
override LDFLAGS += -L.

#------------------------------------------------------------------------------#

CC = gcc
FC = gfortran
FFLAGS ?= -O3 -march=native
# FFLAGS = -g -Og -mieee-fp -fcheck=all
FFLAGS += -fimplicit-none -ffree-line-length-none -pedantic
FFLAGS += -Wall -Wno-unused-dummy-argument -Warray-temporaries -Wrealloc-lhs

# CC = icc
# FC = ifort
# FFLAGS = -O2 -xhost -fp-model precise
# # FFLAGS = -g -Og -mieee-fp -C
# FFLAGS += -warn all,nounused

#------------------------------------------------------------------------------#

# search for sources in these directories
VPATH = ../src:../src/util:../src/prog:../src/math

# objects that go into the shared library
objects = alphadisk.o balance.o fileunits.o globals.o grid.o heatbil.o modelmag.o relaxation_c.o relaxation.o rk4settings.o rxsettings.o ss73solution.o
# various math utilities
objects += bisect.o cgs.o deriv.o eulerintegr.o findzer.o findzer_multi.o histogram.o integrate.o interpol.o linsect.o random.o rk4integr.o threshold.o ramps.o
# objects that are needed for binary programs only
objects_util = results.o settings.o summary.o

#------------------------------------------------------------------------------#

# if using intel compiler, use the MKL instead of standard LAPACK
ifeq ($(FC),ifort)
LDLIBS += -mkl=sequential
else
LDLIBS += -llapack
endif

#------------------------------------------------------------------------------#

programs = dv-alpha dv-alpha-rx dv-mag dv-mag-rx dv-rad1

#------------------------------------------------------------------------------#

all: $(programs)

lib: libdiskvert.so

#------------------------------------------------------------------------------#

install: all
	install -d $(bindir)
	install $(programs) $(bindir)

install-lib: lib
	install -d $(libdir)
	install libdiskvert.so $(libdir)
	install -d $(fmoddir)/diskvert
	install -m 644 slf_{cgs,rk4integr,threshold}.mod \
	 	{globals,settings,fileunits,grid,rk4settings}.mod \
		{alphadisk,modelmag,relaxation,ss73solution}.mod \
		$(fmoddir)/diskvert

install-user: prefix = $(HOME)/.local
install-user: install

install-lib-user: prefix = $(HOME)/.local
install-lib-user: install-lib

#------------------------------------------------------------------------------#

%.o %.mod: %.F90
	$(FC) $(INCLUDE) $(CPPFLAGS) $(FFLAGS) -c $< -o $@
%.o %.mod: %.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c $< -o $@

include deps.inc

relaxation.o    : mrxcoeff.fi mrxdims.fi mrxhash.fi mrxptrs.fi
settings.o      : ../libconfort/libconfort.a

#------------------------------------------------------------------------------#

$(programs): %: $(objects_util) $(objects) %.o ../libconfort/libconfort.a
	$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

#------------------------------------------------------------------------------#

libdiskvert.so: override FFLAGS += -fPIC
libdiskvert.so: $(objects)
	$(FC) $(LDFLAGS) -shared $^ $(LDLIBS) -o $@

#------------------------------------------------------------------------------#

../libconfort/libconfort.a:
	$(MAKE) -C ../libconfort libconfort.a  \
				CC="$(CC)" CFLAGS="-g -O2" \
				FC="$(FC)" FFLAGS="-g -O2"

#------------------------------------------------------------------------------#

clean:
	$(RM) *.mod *.smod *.a *.o $(programs) *.so diskvert-paths.sh
	$(MAKE) -C ../libconfort distclean
