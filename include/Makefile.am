BUILT_SOURCES = .done-all-links

AUTOMAKE_OPTIONS = -Wno-portability

DIRLINKS = $(top_srcdir)/Repository \
           $(top_srcdir)/EventRecord \
           $(top_srcdir)/Pointer \
           $(top_srcdir)/StandardModel \
           $(top_srcdir)/Vectors \
           $(top_srcdir)/PDT \
           $(top_srcdir)/ACDC \
           $(top_srcdir)/Interface \
           $(top_srcdir)/Handlers \
           $(top_srcdir)/LesHouches \
           $(top_srcdir)/Persistency \
           $(top_srcdir)/PDF \
           $(top_srcdir)/Helicity \
           $(top_srcdir)/Utilities \
           $(top_srcdir)/MatrixElement \
           $(top_srcdir)/Cuts \
           $(top_srcdir)/Analysis

CONFIGHEADERS = $(top_srcdir)/Config/Complex.h \
		$(top_srcdir)/Config/Constants.h \
                $(top_srcdir)/Config/Containers.h \
		$(top_srcdir)/Config/Pointers.h \
                $(top_srcdir)/Config/ThePEG.h \
                $(top_srcdir)/Config/algorithm.h \
		$(top_srcdir)/Config/PhysicalQty.h \
		$(top_srcdir)/Config/PhysicalQtyComplex.h \
		$(top_srcdir)/Config/PhysicalQtyOps.h \
		$(top_srcdir)/Config/TemplateTools.h \
		$(top_srcdir)/Config/Unitsystem.h \
		$(top_srcdir)/Config/HepMCHelper.h \
		$(top_srcdir)/Config/std.h

CLEANFILES = .done-all-links

.done-all-links: $(DIRLINKS) $(CONFIGHEADERS)
@EMPTY@ifdef SHOWCOMMAND
	mkdir -p ThePEG/Config
	$(LN_S) -f $(addprefix ../, $(DIRLINKS)) ThePEG
	$(LN_S) -f $(addprefix ../../, $(CONFIGHEADERS)) ThePEG/Config
	touch .done-all-links
@EMPTY@else
	@echo "sym-linking header files..."
	@mkdir -p ThePEG/Config
	@$(LN_S) -f $(addprefix ../, $(DIRLINKS)) ThePEG
	@$(LN_S) -f $(addprefix ../../, $(CONFIGHEADERS)) ThePEG/Config
	@touch .done-all-links
@EMPTY@endif

install-data-local:
	find ThePEG -follow \( -name '*.h' -or -name '*.icc' \
	-or -name '*.tcc' -or -name '*.fh' -or -name '*.xh' \) \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/ThePEG

clean-local:
	rm -rf ThePEG
