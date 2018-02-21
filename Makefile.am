if JAVAGUI
  JAVADIR = java
endif
MAINDIRS = include Config Utilities Helicity Interface LesHouches Vectors \
          PDT PDF Persistency Handlers MatrixElement Pointer ACDC \
          StandardModel Repository EventRecord Cuts Analysis lib src Doc

DIST_SUBDIRS = $(MAINDIRS) java

SUBDIRS = $(MAINDIRS) $(JAVADIR)

EXTRA_DIST = GUIDELINES

DISTCHECK_CONFIGURE_FLAGS = --without-lhapdf --with-gsl=$(GSLPATH)

ACLOCAL_AMFLAGS = -I m4

DISTCLEANFILES = config.thepeg
