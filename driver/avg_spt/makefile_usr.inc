$(OBJDIR)/avg_spt.o        :avg_spt.f;                         $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/avg_spt_usr.o    :avg_spt_usr.f;                     $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/avg_spt_debug.o  :avg_spt_debug.f;                     $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@

