include ../../make.inc
# Stolen from https://web.mit.edu/gnu/doc/html/make_4.html
objects := $(patsubst %.f,%.o,$(wildcard *.f))
# Groups
all: lib $(objects)

clean:
	rm *.o myRefLib.a

lib: myRefLib.a

# Compile all fortran files the same way
$(objects): %.o: %.f
	$(FC) -c $(FFLAGS) -O3 $< -o $@

# Compile all our objects into a single library
myRefLib.a: $(objects)
	ar rcs myRefLib.a $(objects)
