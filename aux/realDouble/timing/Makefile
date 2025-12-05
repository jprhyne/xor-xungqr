include ../make.inc

.PHONY: all clean
all: time_dlarft.exe

clean:
	+$(MAKE) -C referenceFiles clean
	rm *.o *.exe

time_dlarft.exe: time_dlarft.o $(MY_LIBS)/myLapack.a referenceFiles/myRefLib.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^  $(ALL_OPT_LIB) -o $@

$(MY_LIBS)/myLapack.a:
	./compile.sh

referenceFiles/myRefLib.a:
	+$(MAKE) -C referenceFiles all

time_dlarft.o: time_dlarft.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

