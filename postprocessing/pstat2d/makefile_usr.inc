$(OBJDIR)/pstat.o           :pstat.f PSTATD FRAMELP;                 $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/pstat_IO.o        :pstat_IO.f PSTATD FRAMELP;              $(F77) -c $(FL2) -I./ $< -o $@