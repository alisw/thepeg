BUILT_SOURCES = done-all-links

AUTOMAKE_OPTIONS = -Wno-portability

DIRLINKS = $(top_srcdir)/BE \
           $(top_srcdir)/Hadronization \
           $(top_srcdir)/Strategy

CONFIGHEADERS = $(top_srcdir)/Config/TheP8I.h $(top_srcdir)/Config/Pythia8Interface.h  $(top_srcdir)/Config/RndmEngine.h

CLEANFILES = done-all-links

done-all-links: $(DIRLINKS) $(CONFIGHEADERS)
@EMPTY@ifdef SHOWCOMMAND
	mkdir -p TheP8I/Config
	$(LN_S) -f $(addprefix ../, $(DIRLINKS)) TheP8I
	$(LN_S) -f $(addprefix ../../, $(CONFIGHEADERS)) TheP8I/Config
	$(LN_S) -f ../../../Config/config.h TheP8I/Config
	touch done-all-links
@EMPTY@else
	@echo "sym-linking header files..."
	@mkdir -p TheP8I/Config
	@$(LN_S) -f $(addprefix ../, $(DIRLINKS)) TheP8I
	@$(LN_S) -f $(addprefix ../../, $(CONFIGHEADERS)) TheP8I/Config
	@$(LN_S) -f ../../../Config/config.h TheP8I/Config
	@touch done-all-links
@EMPTY@endif

install-data-local:
	find TheP8I -follow \( -name '*.h' -or -name '*.icc' \
	-or -name '*.tcc' -or -name '*.fh' -or -name '*.xh' \) \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/TheP8I
	if test -x $(DESTDIR)$(bindir)/setupThePEG; then $(DESTDIR)$(bindir)/setupThePEG -init -r $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo -o $(DESTDIR)$(pkglibdir)/ThePEGDefaults.rpo $(top_srcdir)/src/TheP8IRemove.in; fi

clean-local:
	rm -rf TheP8I
