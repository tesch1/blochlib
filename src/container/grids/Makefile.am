
noinst_LTLIBRARIES = libblochlibmaind.la


include_HEADERS = coords.h unigrid.h gengrid.h fullcubegrid.h grids.h fullspheregrid.h surfspheregrid.h surfcubegrid.h xyzgenshape.h xyzshapeparts.h xyzshape.h edgedetect.h

libblochlibmaind_la_SOURCES = unigrid.cc fullcubegrid.cc fullspheregrid.cc surfspheregrid.cc surfcubegrid.cc xyzshapeparts.cc xyzshape.cc


INCLUDES = -I../../ @FFTW_INC@

AM_CXXFLAGS = @CXXFLAGS@

includedir=${prefix}/include/blochlib/container/grids/

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

clean:
	rm -rf .libs _libs
	rm -f *.o core *.core
	rm -f *.lo
