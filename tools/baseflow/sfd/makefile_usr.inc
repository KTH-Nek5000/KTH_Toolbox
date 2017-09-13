$(OBJDIR)/sfd.o             :sfd.f SFD CHKPTMSTPD CHKPOINTD;         $(F77) -c $(FL2) -I./ $< -o $@
