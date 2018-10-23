


noinst_LTLIBRARIES = libblochlibmaind.la

#EXTRA_DIST = ham_gen.cc 

include_HEADERS =  powder.h compute.h oneFID.h ham_gen.h 

libblochlibmaind_la_SOURCES = powder.cc ham_gen.cc

INCLUDES = -I../ @FFTW_INC@

#AM_CXXFLAGS = @CXXFLAGS@


includedir=${prefix}/include/blochlib/QMspins/driver/

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



