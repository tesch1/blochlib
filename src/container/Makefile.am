
noinst_LTLIBRARIES = libblochlibmaind.la



EXTRA_DIST =gaussjordon.cc MasterHeader.h 

include_HEADERS = complex.h const_itter.h MemChunk.h operations.h rankType.h vector_add.h containers.h range.h

SUBDIRS = matrix Vector grids 

libblochlibmaind_la_SOURCES = complex.cc range.cc


INCLUDES = -I../

AM_CXXFLAGS = @CXXFLAGS@

includedir=${prefix}/include/blochlib/container/

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
