CLEANFILES = AllInterfaces.h MoreInterfaces.h doxywarnings.log

docdir = $(pkgdatadir)/Doc

doc_DATA = refman.conf refman.h

dist_doc_DATA = repository.hlp MakeDocs.in

dist_noinst_SCRIPTS = fixinterfaces.pl

clean-local: html-clean

html-clean:
	rm -rf refman-html ThePEG-refman.tag

html-local: refman-html

AllInterfaces.h: MakeDocs.in
	$(top_builddir)/src/setupThePEG --init --exitonerror -L$(top_builddir)/lib $(srcdir)/MakeDocs.in

MoreInterfaces.h: MakeDocs.in
	$(top_builddir)/src/setupThePEG --init --exitonerror -L$(top_builddir)/lib $(srcdir)/MakeDocs.in

refman-html: refman.conf refman.h AllInterfaces.h MoreInterfaces.h
	$(builddir)/fixinterfaces.pl AllInterfaces.h
	$(builddir)/fixinterfaces.pl MoreInterfaces.h
	doxygen $<

install-data-local:
	if test -d refman-html; then $(mkdir_p) $(DESTDIR)$(docdir);\
	find refman-html -follow -name '*.*' \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(docdir)/\{\} \;; \
	$(install_sh_DATA) ThePEG-refman.tag $(DESTDIR)$(docdir)/ThePEG-refman.tag; fi

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)/refman-html \
	$(DESTDIR)$(docdir)/ThePEG-refman.tag
