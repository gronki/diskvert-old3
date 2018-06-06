VERSION = 180302

#------------------------------------------------------------------------------#

prefix = /usr/local
bindir = $(prefix)/bin
libdir = $(prefix)/lib
fmoddir = $(libdir)/gfortran/modules

#------------------------------------------------------------------------------#

override INCLUDE += -I../libconfort
override LDFLAGS += -L.

#------------------------------------------------------------------------------#

CC = gcc
FC = gfortran

# release
FFLAGS += -O3 -flto
# debug
# FFLAGS := -g -Og -fcheck=all -fbacktrace -mieee-fp

FFLAGS += -fimplicit-none -ffree-line-length-none -pedantic
FFLAGS += -Wall -Wno-unused-dummy-argument -Warray-temporaries -Wrealloc-lhs

#------------------------------------------------------------------------------#

# CC = icc
# FC = ifort
#
# # release
# FFLAGS := -O3 -xhost -ipo
# # debug
# # FFLAGS := -g -Og -C -mieee-fp
#
# FFLAGS += -warn

#------------------------------------------------------------------------------#

# search for sources in these directories
VPATH = ../src:../src/util:../src/prog:../src/math

# objects that go into the shared library
objects = balance.o fileunits.o globals.o grid.o heatbil.o ss73solution.o
objects += rk4settings.o modelmag.o alphadisk.o
objects += relaxation.o rxsettings.o
# various math utilities
objects += cgs.o threshold.o rk4integr.o eulerintegr.o ramps.o integrate.o interpol.o deriv.o
# objects += bisect.o findzer.o findzer_multi.o histogram.o linsect.o random.o

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

lib: override FFLAGS += -fPIC
lib: libdiskvert.so all

#------------------------------------------------------------------------------#

install: all
	install -d $(DESTDIR)$(bindir)
	install $(programs) $(DESTDIR)$(bindir)

install-lib: lib install
	install -d $(DESTDIR)$(libdir)
	install libdiskvert.so $(DESTDIR)$(libdir)
	install -d $(DESTDIR)$(fmoddir)/diskvert
	install -m 644 slf_{cgs,rk4integr,threshold}.mod \
	 	{globals,settings,fileunits,grid,rk4settings}.mod \
		{alphadisk,modelmag,relaxation,ss73solution}.mod \
		$(DESTDIR)$(fmoddir)/diskvert

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

$(programs): %: $(objects) results.o settings.o %.o ../libconfort/libconfort.a
	$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

#------------------------------------------------------------------------------#

libdiskvert.so: $(objects) relaxation_c.o
	$(FC) $(LDFLAGS) -shared $^ $(LDLIBS) -o $@

#------------------------------------------------------------------------------#

../libconfort/libconfort.a:
	$(MAKE) -C ../libconfort CC=$(CC) FC=$(FC) CFLAGS='-g -O2' FFLAGS='-g -O2' libconfort.a

#------------------------------------------------------------------------------#

clean:
	$(RM) *.mod *.smod *.a *.o
	$(RM) $(programs) *.so
	$(MAKE) -C ../libconfort distclean
