$(OBJDIR)/chkptd.o          :$(CHKPTDMM_SRC)/chkptd.f;               $(F77) -c $(FL2) $< -o $@
