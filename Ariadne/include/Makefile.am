BUILT_SOURCES = done-all-links

AUTOMAKE_OPTIONS = -Wno-portability

DIRLINKS = $(top_srcdir)/DipoleCascade $(top_srcdir)/Cascade $(top_srcdir)/DIPSY

CONFIGHEADERS = $(top_srcdir)/Config/Ariadne.h \
                $(top_srcdir)/Config/Ariadne5.h \
                $(top_srcdir)/Config/CloneBase.fh \
                $(top_srcdir)/Config/CloneBase.h \
                $(top_srcdir)/Config/CloneBase.icc \
                $(top_srcdir)/Config/UnitFO.h

CLEANFILES = done-all-links

done-all-links: $(CONFIGHEADERS) $(DIRLINKS)
@EMPTY@ifdef SHOWCOMMAND
	mkdir -p Ariadne/Config
	$(LN_S) -f $(addprefix ../, $(DIRLINKS)) Ariadne
	$(LN_S) -f $(addprefix ../../, $(CONFIGHEADERS)) Ariadne/Config
	$(LN_S) -f ../../../Config/config.h Ariadne/Config
	touch done-all-links
@EMPTY@else
	@echo "sym-linking header files..."
	@mkdir -p Ariadne/Config
	@$(LN_S) -f $(addprefix ../, $(DIRLINKS)) Ariadne
	@$(LN_S) -f $(addprefix ../../, $(CONFIGHEADERS)) Ariadne/Config
	@$(LN_S) -f ../../../Config/config.h Ariadne/Config
	@touch done-all-links
@EMPTY@endif

install-data-local:
	find Ariadne -follow \( -name '*.h' -or -name '*.icc' \
	-or -name '*.tcc' -or -name '*.fh' -or -name '*.xh' \) \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/Ariadne
	if test -x $(DESTDIR)$(bindir)/setupThePEG; then $(DESTDIR)$(bindir)/setupThePEG -init -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo -o $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo $(top_srcdir)/src/AriadneRemove.in; fi

clean-local:
	rm -rf Ariadne
