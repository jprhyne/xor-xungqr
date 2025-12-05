include ../make.inc

.PHONY: all clean
all: test_dstrk.exe test_dst3rk.exe test_dtrmvoop.exe test_dtrmmoop_lvl2.exe

clean:
	rm *.o *.exe

test_dtrmmoop_lvl2.exe: test_dtrmmoop_lvl2.o
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(MY_LIBS)/myLapack.a $(ALL_OPT_LIB) -o $@

test_dtrmmoop_lvl2.o: test_dtrmmoop_lvl2.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

test_dtrmvoop.exe: test_dtrmvoop.o
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(MY_LIBS)/myLapack.a $(ALL_OPT_LIB) -o $@

test_dtrmvoop.o: test_dtrmvoop.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

test_dstrk.exe: test_dstrk.o
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(MY_LIBS)/myLapack.a $(ALL_OPT_LIB) -o $@

test_dstrk.o: test_dstrk.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

test_dst3rk.exe: test_dst3rk.o
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(MY_LIBS)/myLapack.a $(ALL_OPT_LIB) -o $@

test_dst3rk.o: test_dst3rk.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@
