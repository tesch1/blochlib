


noinst_LTLIBRARIES = libblochmaind.la


include_HEADERS = Isotope.h  bloch.h bloch_basic.h bloch_interac.h bloch_interac_meth.h bloch_int_offset.h mag_fields.h bloch_int_bfield.h bloch_int_relax.h bloch_int_dipdip.h bloch_int_raddamp.h bloch_int_demag.h bloch_int_dimlessdemag.h bloch_int_demagmodulate.h bloch_int_bulks.h bloch_int_diffu.h bloch_int_qmdip.h bloch_int_dimlessdip.h blochParams.h blochParamsBasic.h listblochparsGrad.h listblochparsGrad_meth.h blochsolver.h listblochpars.h listblochpars_meth.h gradientgrid.h pulse.h lyapunov.h scalefunc.h gradfunc.h gradimps.h biot_basic.h biot.h biot_meth.h rotating_grid.h

libblochmaind_la_SOURCES = Isotope.cc blochParamsBasic.cc blochParamsDensity.cc blochParamsParticle.cc gradientgrid.cc pulse.cc gradfunc.cc listblochpars_meth.cc lyapunov_meth.cc scalefunc.cc biot.cc biot_basic.cc bloch_int_demagmodulate.cc


INCLUDES = -I../ @MPI_INC@ @FFTW_INC@



AM_CXXFLAGS = @CXXFLAGS@


includedir=${prefix}/include/blochlib/bloch/

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



