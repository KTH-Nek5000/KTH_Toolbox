$(OBJDIR)/arn_arp.o         :$(ARNARP_SRC)/arn_arp.f ./ARN_ARPD;     $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/arn_arp_io.o      :$(ARNARP_SRC)/arn_arp_io.f ./ARN_ARPD;  $(F77) -c $(FL2) $< -o $@
