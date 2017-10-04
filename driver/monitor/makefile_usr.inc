$(OBJDIR)/mntr_block.o      :mntr_block.f MNTRD;                     $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntr.o            :mntr.f MNTRD;                           $(F77) -c $(FL2) -I./ $< -o $@
