$(OBJDIR)/chkpt_mstp.o      :chkpt_mstp.f CHKPTMSTPD;      $(F77) -c $(FL2) -I./ $< -o $@
