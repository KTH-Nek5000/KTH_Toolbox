$(OBJDIR)/rprm_block.o      :$(RPRM_SRC)/rprm_block.f;               $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/rprm.o            :$(RPRM_SRC)/rprm.f;                     $(FC) -c $(FL2) $< -o $@
