JAVASOURCES = SetupThePEG.java ObjectFrame.java \
              Parameter.java FullSlider.java \
              PopCardPanel.java SwitchOption.java \
              Switch.java RefRef.java Reference.java \
              Interface.java RefVector.java \
              BrowserTree.java ObjectNode.java \
              ObjectRef.java ClassSelector.java \
              ObjectSelector.java StringParameter.java \
              Command.java ParVector.java RunSelector.java \
              RunFrame.java

CLEANFILES = ThePEG.jar thepeg.sh

jardir = $(pkglibdir)
nodist_jar_DATA = ThePEG.jar

dist_noinst_DATA = $(JAVASOURCES) jar-manifest

## .sh needed for OS X case-insensitive file systems
## to avoid confusion with ThePEG directory
noinst_SCRIPTS = thepeg.sh
EXTRA_DIST = thepeg.install

clean-local:
	rm -rf ThePEG

ThePEG:
	mkdir -p ThePEG
	for file in $(JAVASOURCES) jar-manifest; do \
           cd ThePEG; $(LN_S) ../$(srcdir)/$$file $$file; cd ..; done

ThePEG.jar: ThePEG $(JAVASOURCES)
	$(JAVAC) -source 1.4 `for file in $(JAVASOURCES); do echo ThePEG/$$file; done`
	$(JAR) cmf ThePEG/jar-manifest ThePEG.jar ThePEG/*.class

thepeg.sh: thepeg.install Makefile
	sed -e s:@pkglibdir[@]:$(pkglibdir):g \
            -e s:@bindir[@]:$(bindir):g \
            -e s:@java[@]:$(JAVA):g $(srcdir)/thepeg.install > thepeg.sh
	chmod +x thepeg.sh


test: ThePEG.jar
	$(JAVA) -jar ThePEG.jar ../src/setupThePEG -d 1 -L../lib \
                -L../../Pythia7/lib -r ../../Pythia7/lib/Pythia7Defaults.rpo

check-local:
	$(JAVA) -Djava.awt.headless=true -jar ThePEG.jar --classcheck

install-exec-hook: thepeg.sh
	$(install_sh) thepeg.sh $(DESTDIR)$(bindir)/thepeg 

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/thepeg

installcheck-local:
	$(DESTDIR)$(bindir)/thepeg --classcheck
