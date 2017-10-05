$(OBJDIR)/rprm_block.o      :rprm_block.f RPRMD;                     $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/rprm.o            :rprm.f RPRMD MNTRLP;                    $(F77) -c $(FL2) -I./ $< -o $@
