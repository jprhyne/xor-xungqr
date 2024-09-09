include ../make.inc
# Stolen from https://web.mit.edu/gnu/doc/html/make_4.html
objects := $(patsubst %.f,%.o,$(wildcard *.f))
# Groups
all: $(objects)

clean:
	rm *.o

# Compile all fortran files the same way
$(objects): %.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@
