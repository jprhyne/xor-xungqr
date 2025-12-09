include ../make.inc

.PHONY: all clean libs
all: libs time_dlarft.exe time_dlarfb0c2.exe time_dorgqr.exe time_dorgqr_k.exe

clean:
	+$(MAKE) -C referenceFiles clean
	rm *.o *.exe

libs:
	+$(MAKE) -C referenceFiles all

time_dlarft.exe: time_dlarft.o referenceFiles/myRefLib.a myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^  $(ALL_OPT_LIB) -o $@

time_dlarfb0c2.exe: time_dlarfb0c2.o referenceFiles/myRefLib.a myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^  $(ALL_OPT_LIB) -o $@

time_dorgqr.exe: time_dorgqr.o referenceFiles/myRefLib.a myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^  $(ALL_OPT_LIB) -o $@

time_dorgqr_k.exe: time_dorgqr_k.o referenceFiles/myRefLib.a myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^  $(ALL_OPT_LIB) -o $@

myLapack.a: $(MY_LIBS)/myLapack.a
	cp $(MY_LIBS)/myLapack.a ./
	ar d myLapack.a dlarft.o zlarft.o clarft.o slarft.o dorgqr.o sorgqr.o dorglq.o sorglq.o

$(MY_LIBS)/myLapack.a:
	myLapackLib

referenceFiles/myRefLib.a:
	+$(MAKE) -C referenceFiles all

time_dlarft.o: time_dlarft.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

time_dlarfb0c2.o: time_dlarfb0c2.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

time_dorgqr.o: time_dorgqr.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

time_dorgqr_k.o: time_dorgqr_k.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

