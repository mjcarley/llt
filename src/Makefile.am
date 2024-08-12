AM_CPPFLAGS = -I$(top_srcdir) -I$(includedir) -DG_LOG_DOMAIN=\"LLT\"

bin_SCRIPTS = llt-xfoil

BUILT_SOURCES = llt-xfoil

llt-xfoil: llt-xfoil.in

lib_LTLIBRARIES = libllt.la

libllt_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic

libllt_la_SOURCES = \
	llt-private.h \
	llt.c \
	lattice.c \
	solver.c \
	assembly.c \
	wings.c \
	sections.c \
	scripts.c \
	naca.c

include_HEADERS = \
	llt.h

CLEANFILES = $(BUILT_SOURCES)
