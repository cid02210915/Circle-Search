bin_PROGRAMS = circles MASTER_circles
circles_SOURCES = inifile.F90 read_archive_map.F90 circles.F90
MASTER_circles_SOURCES = \
    inifile.F90 \
    read_archive_map.F90 \
    MASTER_circles.F90
CLEANFILES = inifile.mod

AM_FCFLAGS+=$(HEALPIX_CFLAGS)
LDADD=$(HEALPIX_LIBS) $(SHARP_LIBS) $(CFITSIO_LIBS)

if WANT_SINGLE_PRECISION
    LDADD+=-lfftw3f
else
    LDADD+=-lfftw3
    AM_FCFLAGS+=-DDFFT
endif

if WANT_MPI
    AM_FCFLAGS+=-DMPI
endif