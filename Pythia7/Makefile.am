SUBDIRS =  include Config Strategy StringFrag Shower lib Doc src

EXTRA_DIST = GUIDELINES

ACLOCAL_AMFLAGS = -I m4

install-recursive: unregister-from-repo

unregister-from-repo:
	@cd src; $(MAKE) $(AM_MAKEFLAGS) unregister-from-repo
