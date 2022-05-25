$(OBJDIR)/stat.o            :$(STAT_SRC)/stat.f ./STATD;             $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/stat_IO.o         :$(STAT_SRC)/stat_IO.f ./STATD;          $(FC) -c $(FL2) $< -o $@
