include ../make.inc

srcObjs := $(patsubst %.f,%.o,$(wildcard ../src/*.f))
depObjs := $(patsubst %.f,%.o,$(wildcard ../../../2r/realSingle/src/*.f))

# Groups
all: objectFiles driver_slarft.exe

clean:
	rm *.o *.exe

objectFiles: 
	$(MAKE) -C ../src all

driver_slarft.exe: driver_slarft.o test_slarft.o $(srcObjs) $(depObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_slarft.o: driver_slarft.c
	$(CC) -c $(CFLAGS) $^ -o $@

test_slarft.o:   test_slarft.f
	$(FC) -c $(FFLAGS) $^ -o $@

