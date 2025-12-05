include ../make.inc

.PHONY: all clean
all: test_zhtrk.exe test_zlarft_ut.exe test_zht3rk.exe test_ztrmvoop.exe test_ztrmmoop_lvl2.exe

clean:
	rm *.o *.exe

test_ztrmvoop.exe: test_ztrmvoop.o
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(MY_LIBS)/myLapack.a $(ALL_OPT_LIB) -o $@

test_ztrmvoop.o: test_ztrmvoop.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

test_zlarft_ut.exe: test_zlarft_ut.o
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(MY_LIBS)/myLapack.a $(ALL_OPT_LIB) -o $@

test_zlarft_ut.o: test_zlarft_ut.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

test_zhtrk.exe: test_zhtrk.o
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(MY_LIBS)/myLapack.a $(ALL_OPT_LIB) -o $@

test_zhtrk.o: test_zhtrk.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

test_zht3rk.exe: test_zht3rk.o
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(MY_LIBS)/myLapack.a $(ALL_OPT_LIB) -o $@

test_zht3rk.o: test_zht3rk.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

test_ztrmmoop_lvl2.exe: test_ztrmmoop_lvl2.o
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(MY_LIBS)/myLapack.a $(ALL_OPT_LIB) -o $@

test_ztrmmoop_lvl2.o: test_ztrmmoop_lvl2.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@
