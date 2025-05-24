include ../make.inc

srcObjs := $(patsubst %.f,%.o,$(wildcard ../src/*.f))
larftObjs :=  ../../../larft/realDouble/src/my_dlarft_rec.o

# Groups
all: objectFiles driver_dtrtrm.exe driver_trmm.exe driver_dlumm.exe driver_dorgkr.exe driver_dorglk.exe driver_dorgkl.exe driver_dorgrk.exe

clean:
	rm *.o *.exe

objectFiles:
	+$(MAKE) -C ../src all
	+$(MAKE) -C ../../../larft/realDouble/src all

driver_dst3rk.exe: driver_dst3rk.o test_dst3rk.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_dtrtrm.exe:   driver_dtrtrm.o test_dtrtrm.o  $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_trmm.exe:   driver_trmm.o testdtrmmoop.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_dorgkr.exe: driver_dorgkr.o test_dorgkr.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_dorglk.exe: driver_dorglk.o test_dorglk.o ../../../larft/realDouble/src/my_dlarft_rec.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_dorgkl.exe: driver_dorgkl.o test_dorgkl.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_dorgrk.exe: driver_dorgrk.o test_dorgrk.o ../../../larft/realDouble/src/my_dlarft_rec.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_dlumm.exe:   driver_dlumm.o test_dlumm.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_dtrtrm.o:   driver_dtrtrm.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_trmm.o:   driver_trmm.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_dlumm.o:   driver_dlumm.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_dorgkr.o: driver_dorgkr.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_dorglk.o: driver_dorglk.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_dorgkl.o: driver_dorgkl.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_dorgrk.o: driver_dorgrk.c
	$(CC) -c $(CFLAGS) $^ -o $@

test_dtrtrm.o:     testdtrtrm.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_dorgkr.o:   test_dorgkr.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_dorglk.o:   test_dorglk.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_dorgkl.o:   test_dorgkl.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_dorgrk.o:   test_dorgrk.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_dlumm.o:     test_dlumm.f
	$(FC) -c $(FFLAGS) $^ -o $@

testdtrmmoop.o:  testdtrmmoop.f
	$(FC) -c $(FFLAGS) $^ -o $@
