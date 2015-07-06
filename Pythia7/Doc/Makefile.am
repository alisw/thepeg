CLEANFILES = AllInterfaces.h Pythia7-refman.tag tagfileThePEG.tag

docdir = $(pkgdatadir)/Doc

dist_doc_DATA = refman.conf refman.h MakeDocs.in

EXTRA_DATA = tagfileThePEG.tag

AllInterfaces.h: MakeDocs.in
	$(SETUPTHEPEG) -L../lib -L$(THEPEGLIB) --exitonerror --init MakeDocs.in

tagfileThePEG.tag:
	@if test -e $(THEPEGDOC)/ThePEG-refman.tag; then $(LN_S) $(THEPEGDOC)/ThePEG-refman.tag tagfileThePEG.tag; fi

clean-local: html-clean

html-clean:
	rm -rf refman-html refman-latex

html-local: refman-html

pdf-local: refman-html
	if test -d refman-latex; then cd refman-latex; make refman.pdf; fi

refman-html: refman.conf refman.h tagfileThePEG.tag AllInterfaces.h
	$(THEPEGDOC)/fixinterfaces.pl AllInterfaces.h
	doxygen $<

install-data-local:
	if test -d refman-html; then $(mkdir_p) $(DESTDIR)$(docdir);\
	find refman-html -follow -name '*.*' \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(docdir)/\{\} \;; \
	$(install_sh_DATA) Pythia7-refman.tag $(DESTDIR)$(docdir)/Pythia7-refman.tag; fi
	if test -d refman-latex -a -f refman-latex/refman.pdf; then \
	$(install_sh_DATA) refman-latex/refman.pdf $(DESTDIR)$(docdir)/refman.pdf; fi

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)/refman-html \
	$(DESTDIR)$(docdir)/Pythia7-refman.tag \
	$(DESTDIR)$(docdir)/refman.pdf

