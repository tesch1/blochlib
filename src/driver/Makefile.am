

include_HEADERS = ckrk.h ckrk_meth.h bs.h bs_meth.h gear.h gear_meth.h integrate.h stiffbs.h stiffbs_meth.h bdf_gear.h bdf_gear_meth.h basicODE.h


#INCLUDES = -I../

#AM_CXXFLAGS = @CXXFLAGS@


includedir=${prefix}/include/blochlib/driver/

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
