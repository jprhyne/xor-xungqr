include ../make.inc

srcObjs := $(patsubst %.f,%.o,$(wildcard ../src/*.f))

# Groups
all: objectFiles driver_ztvt.exe driver_trmm.exe driver_zungkr.exe driver_lumm.exe driver_zlarft.exe driver_zst3rk.exe

clean:
	rm *.o *.exe

objectFiles: 
	$(MAKE) -C ../src all

driver_zst3rk.exe: driver_zst3rk.o test_zst3rk.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_ztvt.exe:   driver_ztvt.o test_ztvt.o  $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_trmm.exe:   driver_trmm.o testztrmmoop.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_zungkr.exe: driver_zungkr.o test_zungkr.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_zlarft.exe: driver_zlarft.o test_zlarft.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_lumm.exe:   driver_lumm.o test_lumm.o $(srcObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_ztvt.o:   driver_ztvt.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_zlarft.o: driver_zlarft.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_trmm.o:   driver_trmm.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_lumm.o:   driver_lumm.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_zungkr.o: driver_zungkr.c
	$(CC) -c $(CFLAGS) $^ -o $@

test_ztvt.o:     test_ztvt.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_zungkr.o:   test_zungkr.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_zlarft.o:   test_zlarft.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_lumm.o:     test_lumm.f
	$(FC) -c $(FFLAGS) $^ -o $@

testztrmmoop.o:  testztrmmoop.f
	$(FC) -c $(FFLAGS) $^ -o $@
