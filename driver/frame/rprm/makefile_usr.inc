$(OBJDIR)/rprm_block.o      :$(RPRM_SRC)/rprm_block.f;               $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/rprm.o            :$(RPRM_SRC)/rprm.f;                     $(F77) -c $(FL2) $< -o $@
