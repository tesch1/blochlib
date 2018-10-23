
noinst_LTLIBRARIES = libblochlibmaind.la

include_HEADERS =  dmatrix.h fmatrix.h genmatrix.h _matrix.h smatrix.h trimatrix.h imatrix.h hmatrix.h  matrix.h matvec.h matmat.h matdiagonalize.h matrixconfig.h matassign.h LU.h matqr.h matexp.h matmatmul.h matfft.h dgemm.h mdmd_a1bc_md.h mdmd_acbc_md.h mdmdt_a1bc_md.h mdmdt_acbc_md.h mdtmd_a1bc_md.h mdtmd_acbc_md.h mdtmdt_a1bc_md.h mdtmdt_acbc_md.h mdmd_cm.h

libblochlibmaind_la_SOURCES = matexp.cc matdiagonalize.cc mdmd_a1bc_md_cm.cc mdmd_acbc_md_cm.cc mdmdt_a1bc_md_cm.cc mdmdt_acbc_md_cm.cc mdtmd_a1bc_md_cm.cc mdtmd_acbc_md_cm.cc mdtmdt_a1bc_md_cm.cc mdtmdt_acbc_md_cm.cc matmatmul.cc


INCLUDES = -I../../ @FFTW_INC@


#AM_CXXFLAGS = @CXXFLAGS@ 

includedir=${prefix}/include/blochlib/container/matrix/

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
