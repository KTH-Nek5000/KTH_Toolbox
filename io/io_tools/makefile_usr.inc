$(OBJDIR)/io_tools_block.o  :io_tools_block.f;   $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_tools.o        :io_tools.f;         $(F77) -c $(FL2) -I./ $< -o $@
