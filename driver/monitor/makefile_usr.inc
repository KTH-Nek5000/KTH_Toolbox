$(OBJDIR)/mntrlog_block.o   :mntrlog_block.f MNTRLOGD;               $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrlog.o         :mntrlog.f MNTRLOGD;                     $(F77) -c $(FL2) -I./ $< -o $@
