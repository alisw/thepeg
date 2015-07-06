CLEANFILES = done-all-links

pkglib_LTLIBRARIES = LEPEventShapes.la MultiplicityCount.la

# Version info should be updated if any interface or persistent I/O
# function is changed

LEPEventShapes_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO)
LEPEventShapes_la_SOURCES = EventShapes.cc EventShapes.fh EventShapes.h\
			 LEPEventShapes.cc LEPEventShapes.fh\
			 LEPEventShapes.h LEPEventShapes.icc

MultiplicityCount_la_LDFLAGS = -module $(LIBTOOLVERSIONINFO)
MultiplicityCount_la_SOURCES = MultiplicityCount.cc MultiplicityCount.fh\
			  MultiplicityCount.h MultiplicityCount.icc

all-local: done-all-links

done-all-links:
@EMPTY@ifdef SHOWCOMMAND
	for file in $(INPUTFILES); do \
	  if test ! -f $$file; then $(LN_S) $(srcdir)/$$file $$file; fi; done
	echo "stamp" > done-all-links
@EMPTY@else
	@echo "sym-linking input files files..."
	@for file in $(INPUTFILES); do \
	  if test ! -f $$file; then $(LN_S) $(srcdir)/$$file $$file; fi; done
	@echo "stamp" > done-all-links
@EMPTY@endif

