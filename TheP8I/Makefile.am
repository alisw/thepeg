if PYTHIA8WORKS

WORKINGDIRS = include Config Strategy Hadronization BE lib Doc src

install-recursive: unregister-from-repo

unregister-from-repo:
	@cd src; $(MAKE) $(AM_MAKEFLAGS) unregister-from-repo

endif

SUBDIRS = $(WORKINGDIRS)

EXTRA_DIST = GUIDELINES

ACLOCAL_AMFLAGS = -I m4

