$(OBJDIR)/tsrs.o              :$(TSRS_SRC)/tsrs.f ./TSRSD;              $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/tsrs_IO.o           :$(TSRS_SRC)/tsrs_IO.f ./TSRSD;           $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/pts_redistribute.o  :$(TSRS_SRC)/pts_redistribute.f ./TSRSD;  $(F77) -c $(FL2) $< -o $@
