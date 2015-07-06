CLEANFILES = AllInterfaces.h

docdir = $(pkgdatadir)/Doc

doc_DATA = refman.conf refman.h

clean-local: html-clean

html-clean:
	rm -rf refman-html refman-latex

html-local: refman-html

pdf-local: refman-html
	if test -d refman-latex; then cd refman-latex; make refman.pdf; fi

refman-html: refman.conf refman.h
	$(THEPEGDOC)/fixinterfaces.pl AllInterfaces.h
	$(LN_S) -f $(THEPEGDOC)/refman.conf ThePEG.conf
	$(LN_S) -f $(THEPEGDOC)/ThePEG-refman.tag .
	$(LN_S) -f $(THEPEGPATH)/include include
	doxygen $<

install-data-local:
	if test -d refman-html; then $(mkdir_p) $(DESTDIR)$(docdir);\
	find refman-html -follow -name '*.*' \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(docdir)/\{\} \;; \
	$(install_sh_DATA) TheP8I-refman.tag $(DESTDIR)$(docdir)/TheP8I-refman.tag; fi
	if test -d refman-latex -a -f refman-latex/refman.pdf; then \
	$(install_sh_DATA) refman-latex/refman.pdf $(DESTDIR)$(docdir)/refman.pdf; fi

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)/refman-html \
	$(DESTDIR)$(docdir)/TheP8I-refman.tag \
	$(DESTDIR)$(docdir)/refman.pdf

