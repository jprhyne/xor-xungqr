include ../make.inc
# Group definitions

objects := $(patsubst %.f,%.o,$(wildcard *.f))

default: $(objects)

clean:
	rm *o 

# Object file creations

$(objects): %.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@
