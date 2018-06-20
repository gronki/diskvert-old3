#------------------------------------------------------------------------------#
#----------------------------- INSTALLATION PATHS -----------------------------#
#------------------------------------------------------------------------------#

prefix = /usr/local
bindir = $(prefix)/bin
libdir = $(prefix)/lib
includedir = $(prefix)/include
fmoddir = $(includedir)

#------------------------------------------------------------------------------#

INCLUDE := -I../libconfort
LDFLAGS := -L.

#------------------------------------------------------------------------------#
#---------------------------- CFLAGS AND COMPILERS ----------------------------#
#------------------------------------------------------------------------------#

# TO USE COMPILERS OTHER THAN GCC, COMMENT THIS BLOCK, AND UNCOMMENT
# ONE OF THE BELOW BLOCKS (Intel or PGI)

CC := gcc
FC := gfortran

# TO SWITCH BETWEEN PERFORMANCE AND DEBUG COMPILATION, UNCOMMENT ONE:
# default
FFLAGS := -O2
# debug (slow)
# FFLAGS := -ggdb -O0 -fcheck=all -mieee-fp

FFLAGS += -fimplicit-none -ffree-line-length-none -fbacktrace
FFLAGS += -Wall -Wno-unused-dummy-argument -pedantic
FFLAGS += -Warray-temporaries -Wrealloc-lhs-all

LDLIBS += -lopenblas

#------------------------------------------------------------------------------#
# Intel compilers

# CC := icc
# FC := ifort
#
# FFLAGS := -O2 -xhost -warn
# LDLIBS += -mkl=sequential

#------------------------------------------------------------------------------#
# PGI compilers (not working)

# CC := pgcc
# FC := pgf90
#
# FFLAGS := -O2
# LDLIBS += -llapack -lblas

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

programs = dv-alpha dv-alpha-rx dv-mag dv-mag-rx dv-rad1

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

all: $(programs)

lib: override FFLAGS += -fPIC
lib: libdiskvert.so

#------------------------------------------------------------------------------#
#------------------------------ INSTALLATION ----------------------------------#
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
#-------------------------------- BUILD RULES ---------------------------------#
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
	$(MAKE) -C ../libconfort CC=$(CC) CFLAGS='-O2' \
		FC=$(FC) FFLAGS='-O2' libconfort.a

#------------------------------------------------------------------------------#
#---------------------------------- CLEANUP -----------------------------------#
#------------------------------------------------------------------------------#

clean:
	$(RM) *.mod *.smod *.a *.o
	$(RM) $(programs) *.so
	$(MAKE) -C ../libconfort distclean

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
