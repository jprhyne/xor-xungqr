include ../make.inc
# Groups
srcObjs := $(patsubst %.f,%.o,$(wildcard ../src/*.f))
refObjs := $(patsubst %.f,%.o,$(wildcard ../referenceFiles/*.f))

default: srcFiles refFiles test_dorgqr.exe

srcFiles: 
	+$(MAKE) -C ../src

refFiles:
	+$(MAKE) -C ../referenceFiles

clean:
	rm *o *.exe

print: 
	echo $(CFLAGS)

# Executables
test_dorgqr.exe: test_dorgqr.o ../referenceFiles/dorgqr_ref.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

# Object files
test_dorgqr.o: test_dorgqr.c
	$(CC) -c $(CFLAGS) $^ -o $@
