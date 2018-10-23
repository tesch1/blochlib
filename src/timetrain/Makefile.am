
noinst_LTLIBRARIES = libblochlibmaind.la



include_HEADERS =  timetrain.h  timetrain_meth.h gentrain.h unitrain.h consttrain.h consttrain_meth.h filetrain.h 

libblochlibmaind_la_SOURCES = gentrain.cc unitrain.cc filetrain.cc


INCLUDES = -I../  @FFTW_INC@
 
#AM_CXXFLAGS = @CXXFLAGS@

includedir=${prefix}/include/blochlib/timetrain/

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
