$(OBJDIR)/mntrlog_block.o   :$(MNTR_SRC)/mntrlog_block.f;            $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/mntrlog.o         :$(MNTR_SRC)/mntrlog.f;                  $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/mntrtmr_block.o   :$(MNTR_SRC)/mntrtmr_block.f;            $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/mntrtmr.o         :$(MNTR_SRC)/mntrtmr.f;                  $(FC) -c $(FL2) $< -o $@
