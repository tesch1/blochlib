
noinst_LTLIBRARIES = libblochlibmaind.la


include_HEADERS = constants.h utils.h blassert.h plotter.h random.h params.h endians.h matlab5.h vnmrstream.h xwinnmr.h spinsight.h parser.h scriptparse.h wavestream.h

libblochlibmaind_la_SOURCES =  constants.cc utils.cc blassert.cc plotter.cc random.cc params.cc endians.cc matlab5.cc vnmrstream.cc xwinnmr.cc spinsight.cc parser.cc scriptparse.cc wavestream.cc

INCLUDES = -I../ @FFTW_INC@

#AM_CXXFLAGS = @CXXFLAGS@

includedir=${prefix}/include/blochlib/utils/

install-includeHEADERS: $(include_HEADERS)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(includedir)
	@list='$(include_HEADERS)'; for p in $$list; do \
	  if test -f "$$p"; then d= ; else d="$(srcdir)/"; fi; \
	  echo " $(INSTALL_DATA) $$d$$p $(DESTDIR)$(includedir)/$$p"; \
	  $(INSTALL_DATA) $$d$$p $(DESTDIR)$(includedir)/$$p; \
	done

uninstall-includeHEADERS:
	@$(NORMAL_UNINSTALL)
	list='$(include_HEADERS)'; for p in $$list; do \
	  rm -f $(DESTDIR)$(includedir)/$$p; \
	done

