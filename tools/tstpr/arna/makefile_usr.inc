$(OBJDIR)/arna.o            :$(ARNA_SRC)/arna.f ./ARNAD ./TSTPRD;    $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/arna_io.o         :$(ARNA_SRC)/arna_io.f ./ARNAD ./TSTPRD; $(F77) -c $(FL2) $< -o $@
