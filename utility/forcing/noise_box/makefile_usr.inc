$(OBJDIR)/noisebx.o         :noisebx.f NIOSEBXD;           $(F77) -c $(FL2) -I./ $< -o $@
