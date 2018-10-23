
AUTOMAKE_OPTIONS = no-dependencies foreign

SUBDIRS = src

EXTRA_DIST = acinclude.m4 bootstrap.sh doc data blcc.in ltmain.sh libtool mkinstalldirs ltconfig WINDOWS.README README.BLAS README.OSX\
crystals matlabfunc benchmark \
examples/README examples/blochsims examples/classes examples/data_readers examples/diffusion examples/mpiplay examples/qmsims examples/windows examples/solid-2.0

AUTOHEADER=echo

#AM_CXXFLAGS =@CXXFLAGS@

AUXSCRIPTS=blcc
AUXBINDIR="$(prefix)/bin/"
install: install-recursive
	chmod uga+x blcc
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(AUXBINDIR)
	@list='$(AUXSCRIPTS)'; for p in $$list; do \
	  if test -f "$$p"; then d= ; else d="$(srcdir)/"; fi; \
	  echo " $(INSTALL_DATA) $$d$$p $(AUXBINDIR)/$$p"; \
	  $(INSTALL_DATA) $$d$$p $(AUXBINDIR)/$$p; \
	  chmod uga+x $(AUXBINDIR)$(AUXSCRIPTS); \
	done



