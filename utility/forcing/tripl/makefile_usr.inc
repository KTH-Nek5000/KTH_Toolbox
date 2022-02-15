$(OBJDIR)/tripl.o           :$(TRIPL_SRC)/tripl.f ./TRIPLD;      $(F77) -c $(FL2) $< -o $@
