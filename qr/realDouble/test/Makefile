include ../make.inc
# Groups
srcObjs := $(patsubst %.f,%.o,$(wildcard ../src/*.f))
refObjs := $(patsubst %.f,%.o,$(wildcard ../referenceFiles/*.f))

default: srcFiles refFiles krFiles test_dorgz.exe

srcFiles: 
	+$(MAKE) -C ../src

refFiles:
	+$(MAKE) -C ../referenceFiles

krFiles:
	+$(MAKE) -C ../../../2r/realDouble/src

clean:
	rm *o *.exe

print: 
	echo $(CFLAGS)

# Executables
test_dorgz.exe: test_dorgz.o myNorm.o compute_errors.o perf_calc.o print_helpers.o $(refObjs) qToL.o qToR.o ../src/dlarfb0c2.o ../../../2r/realDouble/src/dtrmmoop.o  ../../../2r/realDouble/src/dtvt.o ../../../2r/realDouble/src/lumm.o ../../../2r/realDouble/src/my_dorgkr.o ../src/dorgqr_dlarfb0c2.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

# Object files
test_dorgz.o: test_dorgz.c
	$(CC) -c $(CFLAGS) $^ -o $@

myNorm.o: myNorm.c
	$(CC) -c $(CFLAGS) $^ -o $@

compute_errors.o: compute_errors.c
	$(CC) -c $(CFLAGS) $^ -o $@

perf_calc.o: perf_calc.c
	$(CC) -c $(CFLAGS) $^ -o $@

print_helpers.o: print_helpers.c
	$(CC) -c $(CFLAGS) $^ -o $@

qToL.o: qToL.f
	$(FC) -c $(FFLAGS) $^ -o $@

qToR.o: qToR.f
	$(FC) -c $(FFLAGS) $^ -o $@
