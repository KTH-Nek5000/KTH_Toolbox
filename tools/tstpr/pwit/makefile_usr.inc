$(OBJDIR)/pwit.o            :$(PWIT_SRC)/pwit.f ./TSTPRD;            $(F77) -c $(FL2) $< -o $@
