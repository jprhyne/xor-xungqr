include ../make.inc
# Groups
srcObjs := $(patsubst %.f,%.o,$(wildcard ../src/*.f))
depObjs := $(patsubst %.f,%.o,$(wildcard ../../../2r/realDouble/src/*.f))

all: objectFiles timeDlarft.exe timeDgeqrf.exe

clean:
	rm *.o *.exe

objectFiles:
	$(MAKE) -C ../src all

timeDlarft.exe: timeDlarft.o $(srcObjs) $(depObjs) ../referenceFiles/dlarft.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

timeDgeqrf.exe: timeDgeqrf.o $(srcObjs) $(depObjs) ../referenceFiles/dgeqrf_rec.o ../referenceFiles/dgeqrf_ref.o ../referenceFiles/dlarft.o ../referenceFiles/dgeqrf_qr3.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

timeDlarft.o: timeDlarft.c
	$(CC) -c $(CFLAGS) $^ -o $@

timeDgeqrf.o: timeDgeqrf.c
	$(CC) -c $(CFLAGS) $^ -o $@
