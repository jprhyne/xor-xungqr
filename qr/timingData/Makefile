include make.inc

refObjs := $(patsubst %.f,%.o,$(wildcard ./referenceFiles/*.f))
objects := $(patsubst %.f,%.o,$(wildcard ./src/*.f))
# Group definitions
default: refFiles main.exe

refFiles:
	+$(MAKE) -C referenceFiles

clean: cleanRef
	rm *.exe
	rm *.o

cleanRef:
	+$(MAKE) -C referenceFiles clean

# Executable compilations
main.exe: main.o $(refObjs) $(objects)
	g++ $(CFLAGS) $^ $(OPTBLAS) -o $@

# Object file creations

main.o: main.c
	$(CC) -c $(CFLAGS) $^ -o $@
