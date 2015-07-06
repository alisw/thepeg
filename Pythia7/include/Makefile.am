BUILT_SOURCES = done-all-links

AUTOMAKE_OPTIONS = -Wno-portability

DIRLINKS = $(top_srcdir)/Strategy \
           $(top_srcdir)/StringFrag \
           $(top_srcdir)/Shower

CONFIGHEADERS = $(top_srcdir)/Config/Pythia7.h

CLEANFILES = done-all-links

done-all-links: $(DIRLINKS) $(CONFIGHEADERS)
@EMPTY@ifdef SHOWCOMMAND
	mkdir -p Pythia7/Config
	$(LN_S) -f $(addprefix ../, $(DIRLINKS)) Pythia7
	$(LN_S) -f $(addprefix ../../, $(CONFIGHEADERS)) Pythia7/Config
	$(LN_S) -f ../../../Config/config.h Pythia7/Config
	touch done-all-links
@EMPTY@else
	@echo "sym-linking header files..."
	@mkdir -p Pythia7/Config
	@$(LN_S) -f $(addprefix ../, $(DIRLINKS)) Pythia7
	@$(LN_S) -f $(addprefix ../../, $(CONFIGHEADERS)) Pythia7/Config
	@$(LN_S) -f ../../../Config/config.h Pythia7/Config
	@touch done-all-links
@EMPTY@endif

install-data-local:
	find Pythia7 -follow \( -name '*.h' -or -name '*.icc' \
	-or -name '*.tcc' -or -name '*.fh' -or -name '*.xh' \) \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/Pythia7
	if test -x $(DESTDIR)$(bindir)/setupThePEG; then $(DESTDIR)$(bindir)/setupThePEG -init -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo -o $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo $(top_srcdir)/src/Pythia7Remove.in; fi

clean-local:
	rm -rf Pythia7
