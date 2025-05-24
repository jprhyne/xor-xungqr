include ../make.inc
# Groups
srcObjs := $(patsubst %.f,%.o,$(wildcard ../src/*.f))
refObjs := $(patsubst %.f,%.o,$(wildcard ../referenceFiles/*.f))
0c2Objs := $(patsubst %.f,%.o,$(wildcard ../src/*0c2.f))
helperObjs := myNorm.o compute_errors.o perf_calc.o print_helpers.o qToL.o qToR.o
#krObjs := ../../../2r/realDouble/src/dtrmmoop.o  ../../../2r/realDouble/src/dtvt.o ../../../2r/realDouble/src/lumm.o ../../../2r/realDouble/src/my_zungkr.o

default: srcFiles refFiles test_zungz.exe #krFiles 

srcFiles: 
	+$(MAKE) -C ../src

refFiles:
	+$(MAKE) -C ../referenceFiles

#krFiles:
#	+$(MAKE) -C ../../../2r/realDouble/src

clean:
	rm *o *.exe

print: 
	echo $(CFLAGS)

# Executables
test_zungz.exe: test_zungz.o $(0c2Objs) $(helperObjs) $(refObjs) # $(krObjs) 
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

# Object files
test_zungz.o: test_zungz.c
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
