


#noinst_LTLIBRARIES = libblochlibmaind.la


include_HEADERS = stencilprep.h stencilbasic.h stencilextent.h stencilmap.h stencilexpr.h stencilvector.h stencilprep_func.h boundary_condition.h boundary_general.h boundary_dirichlet.h boundary_neumann.h boundary_extrap.h boundary_edgeblend.h boundary_cornerblend.h


#libblochlibmaind_la_SOURCES = 


INCLUDES = -I../ @FFTW_INC@

#AM_CXXFLAGS = @CXXFLAGS@


includedir=${prefix}/include/blochlib/stencils/

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



