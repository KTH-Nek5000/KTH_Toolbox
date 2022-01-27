$(OBJDIR)/gSyEM.o           :$(GSYEM_SRC)/gSyEM.f ./GSYEMD;          $(F77) -c $(FL2) $< -o $@
