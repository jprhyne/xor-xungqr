include ../make.inc

srcObjs := $(patsubst %.f,%.o,$(wildcard ../src/*.f))
depObjs := $(patsubst %.f,%.o,$(wildcard ../../../2r/realDouble/src/*.f))

# Groups
all: objectFiles driver_dlarft.exe

clean:
	rm *.o *.exe

objectFiles: 
	$(MAKE) -C ../src all
	$(MAKE) -C ../referenceFiles all

driver_dlarft.exe: driver_dlarft.o test_dlarft.o ../referenceFiles/dlarft.o $(srcObjs) $(depObjs)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver_dlarft.o: driver_dlarft.c
	$(CC) -c $(CFLAGS) $^ -o $@

test_dlarft.o:   test_dlarft.f
	$(FC) -c $(FFLAGS) $^ -o $@

