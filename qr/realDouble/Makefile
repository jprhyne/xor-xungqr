include make.inc

2rObjs := $(patsubst %.f,%.o,$(wildcard ../org2r/src/*.f))
objects := $(patsubst %.f,%.o,$(wildcard ./src/*.f))
# Group definitions
default: srcFiles test.exe

srcFiles:
	+$(MAKE) -C src

clean: cleanSrc
	rm driver *exe

cleanSrc:
	+$(MAKE) -C src clean

# Executable compilations
test.exe: test.o $(2rObjs) $(objects)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

# Object file creations

test.o: test.c
	$(CC) -c $(CFLAGS) $^ -o $@
