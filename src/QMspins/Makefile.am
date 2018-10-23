


noinst_LTLIBRARIES = libblochlibmaind.la


include_HEADERS = multispinop.h singlespinop.h spinsys.h spinsyspars.h spin_ten.h space_ten.h j.h dip.h csa.h qua.h sys.h

libblochlibmaind_la_SOURCES = multispinop.cc singlespinop.cc spinsys.cc spinsyspars.cc spin_ten.cc space_ten.cc j.cc qua.cc csa.cc dip.cc sys.cc

SUBDIRS = driver

INCLUDES = -I../ @FFTW_INC@

#AM_CXXFLAGS = @CXXFLAGS@


includedir=${prefix}/include/blochlib/QMspins/

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



