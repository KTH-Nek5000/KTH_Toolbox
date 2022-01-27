$(OBJDIR)/trip.o            :$(TRIP_LINE_SRC)/trip.f ./TRIPD;        $(F77) -c $(FL2) $< -o $@
