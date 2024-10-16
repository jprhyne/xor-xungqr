include ../make.inc
# Groups

srcObjs := $(patsubst %.f,%.o,$(wildcard ../src/*.f))
depObjs := $(patsubst %.f,%.o,$(wildcard ../../../2r/realSingle/src/*.f))

all: objectFiles timeSlarft.exe

clean:
	rm *.o *.exe

objectFiles:
	$(MAKE) -C ../src all

timeSlarft.exe: timeSlarft.o $(srcObjs) $(depObjs) ../referenceFiles/slarft.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

timeSlarft.o: timeSlarft.c
	$(CC) -c $(CFLAGS) $^ -o $@
