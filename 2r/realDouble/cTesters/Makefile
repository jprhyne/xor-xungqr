include ../make.inc
# Groups

srcObjs := $(patsubst %.f,%.o,$(wildcard ../src/*.f))
larftObjs := $(patsubst %.f,%.o,$(wildcard ../../larft/realDouble/src/*.f))

all: objectFiles larftObjects overallTest_optBlas.exe overallTest.exe test_dtrmmoop.exe test_dtrmmoop_optBlas.exe timeGemm.exe timeDlarft.exe

clean:
	rm *.o *.exe

objectFiles:
	$(MAKE) -C ../src all

objectFiles:
	$(MAKE) -C ../../larft/realDouble/src all

overallTest_optBlas.exe: overallTest_opt.o $(srcObjs) ../referenceFiles/dorg2r.o
	$(FC) -D USE_AOCL $(CFLAGS) $^ $(OPTBLAS) -o $@

overallTest.exe: overallTest.o $(srcObjs) ../referenceFiles/dorg2r.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(CBLASLIB) $(BLASLIB) -o $@

test_dtrmmoop.exe: test_dtrmmoop.o  ../src/dtrmmoop.o
	$(FC) $^ $(LAPACKLIB) $(BLASLIB) -o $@

test_dtrmmoop_optBlas.exe: test_dtrmmoop_opt.o ../src/dtrmmoop.o
	$(FC) $^ $(OPTBLAS) -o $@

overallTest.o: overallTest.c 
	$(CC) -c $(CFLAGS) $^ -o $@

overallTest_opt.o: overallTest.c 
	$(CC) -D USE_AOCL -c $(CFLAGS) -I$(OPTBLASINC) -I$(OPTLAPACKINC) $^ -o $@

test_dtrmmoop.o: test_dtrmmoop.c
	$(FC) -c $(CFLAGS) $^ -o $@

test_dtrmmoop_opt.o: test_dtrmmoop.c
	$(CC) -D USE_AOCL -c $(CFLAGS) $^ -o $@

timeGemm.exe: timeGemm.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

timeGemm.o: timeGemm.c
	$(CC) -c $(CFLAGS) $^ -o $@

timeDlarft.exe: timeDlarft.o $(srcObjs) ../referenceFiles/dorg2r.o ../referenceFiles/dlarft.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

timeDlarft.o: timeDlarft.c
	$(CC) -c $(CFLAGS) $^ -o $@
