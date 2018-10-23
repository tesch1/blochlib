
noinst_LTLIBRARIES = libblochlibmaind.la

include_HEADERS = mpi_packer.h mpi_config.h mpi_tools.h mpi_controller.h

libblochlibmaind_la_SOURCES = mpi_config.cc mpi_tools.cc mpi_controller.cc


INCLUDES = -I../ @MPI_INC@ @FFTW_INC@

#AM_CXXFLAGS = @CXXFLAGS@

includedir=${prefix}/include/blochlib/mpi/

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
