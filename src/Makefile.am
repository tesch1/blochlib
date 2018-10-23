lib_LTLIBRARIES = libbloch.la


include_HEADERS = blochlib.h blochconfig.h version.h

EXTRA_DIST = blochconfig.h.in blochconfigCW.h

SUBDIRS = mpi minuit driver container QMspins stencils utils bloch timetrain 

libbloch_la_SOURCES = version.cc  

libbloch_la_LIBADD =  container/complex.lo container/range.lo utils/blassert.lo utils/params.lo utils/endians.lo utils/matlab5.lo utils/vnmrstream.lo utils/xwinnmr.lo utils/spinsight.lo utils/utils.lo utils/constants.lo utils/random.lo utils/plotter.lo utils/parser.lo utils/scriptparse.lo utils/wavestream.lo mpi/mpi_config.lo mpi/mpi_tools.lo mpi/mpi_controller.lo container/matrix/matmatmul.lo container/matrix/matdiagonalize.lo container/matrix/matexp.lo container/grids/fullcubegrid.lo container/grids/unigrid.lo container/grids/surfspheregrid.lo container/grids/fullspheregrid.lo container/grids/surfcubegrid.lo container/grids/xyzshape.lo container/grids/xyzshapeparts.lo  timetrain/filetrain.lo timetrain/gentrain.lo timetrain/unitrain.lo bloch/Isotope.lo bloch/blochParamsBasic.lo bloch/blochParamsParticle.lo bloch/blochParamsDensity.lo bloch/pulse.lo bloch/gradientgrid.lo bloch/listblochpars_meth.lo bloch/bloch_int_demagmodulate.lo bloch/lyapunov_meth.lo bloch/gradfunc.lo bloch/scalefunc.lo bloch/biot_basic.lo bloch/biot.lo QMspins/multispinop.lo QMspins/singlespinop.lo QMspins/spinsys.lo QMspins/spinsyspars.lo QMspins/space_ten.lo QMspins/spin_ten.lo QMspins/sys.lo QMspins/j.lo QMspins/csa.lo QMspins/dip.lo QMspins/qua.lo QMspins/driver/powder.lo QMspins/driver/ham_gen.lo @MINUIT_LIB@ @BLASLIB@

#container/matrix/mdmd_a1bc_md_cm.lo container/matrix/mdmd_acbc_md_cm.lo container/matrix/mdmdt_a1bc_md_cm.lo container/matrix/mdmdt_acbc_md_cm.lo container/matrix/mdtmd_a1bc_md_cm.lo container/matrix/mdtmd_acbc_md_cm.lo container/matrix/mdtmdt_a1bc_md_cm.lo container/matrix/mdtmdt_acbc_md_cm.lo  

#AM_CXXFLAGS = @CXXFLAGS@

AUTOHEADER=echo

libbloch_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

includedir=${prefix}/include/blochlib/


install-include_HEADERS: $(include_HEADERS)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(includedir)
	@list='$(include_HEADERS)'; for p in $$list; do \
	  if test -f "$$p"; then d= ; else d="$(srcdir)/"; fi; \
	  echo " $(INSTALL_DATA) $$d$$p $(DESTDIR)$(includedir)/$$p"; \
	  $(INSTALL_DATA) $$d$$p $(DESTDIR)$(includedir)/$$p; \
	done

uninstall-include_HEADERS:
	@$(NORMAL_UNINSTALL)
	list='$(include_HEADERS)'; for p in $$list; do \
	  rm -f $(DESTDIR)$(includedir)/$$p; \
	done

	
#all:
#	rm -f libblochlib.a
#	ar r libblochlib.a $(OTHEROBJS)
	
	
	

	
