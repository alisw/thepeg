LINKDIRS = Strategy StringFrag Shower

all-local: done-all-links

done-all-links:
@EMPTY@ifdef SHOWCOMMAND
	for dir in $(LINKDIRS); do \
	  for sofile in ../$$dir/.libs/*.so* ../$$dir/.libs/*.la; \
            do if test -f $$sofile; then $(LN_S) -f $$sofile .; fi; done; done
	echo "stamp" > done-all-links
@EMPTY@else
	@echo "Sym-linking all shared libraries ..."
	@for dir in $(LINKDIRS); do \
	  for sofile in ../$$dir/.libs/*.so* ../$$dir/.libs/*.la; \
            do if test -f $$sofile; then $(LN_S) -f $$sofile .; fi; done; done
	@echo "stamp" > done-all-links
@EMPTY@endif

clean-local:
	rm -f *.so* *.la *.rpo done-all-links

