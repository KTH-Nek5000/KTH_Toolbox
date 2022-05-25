$(OBJDIR)/pstat2D.o         :$(PSTAT2D_SRC)/pstat2D.f;               $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/pstat2D_IO.o      :$(PSTAT2D_SRC)/pstat2D_IO.f;            $(FC) -c $(FL2) $< -o $@
