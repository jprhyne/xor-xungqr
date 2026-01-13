include ../make.inc
BINDIR=bin
OBJDIR=obj
SRCDIR=src

.PHONY: all clean libs tex
all: libs $(BINDIR)/count_trtrm.exe $(BINDIR)/count_lumm.exe

clean:
	+$(MAKE) -C referenceFiles clean
	rm $(OBJDIR)/*.o $(BINDIR)/*.exe

libs: myLapack.a
	+$(MAKE) -C referenceFiles all

$(BINDIR)/count_trtrm.exe: $(OBJDIR)/count_trtrm.o referenceFiles/myRefLib.a myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(ALL_OPT_LIB) -o $@

$(BINDIR)/count_lumm.exe: $(OBJDIR)/count_lumm.o referenceFiles/myRefLib.a myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(ALL_OPT_LIB) -o $@

referenceFiles/myRefLib.a:
	+$(MAKE) -C referenceFiles all

myLapack.a: $(MY_LIBS)/myLapack.a
	cp $(MY_LIBS)/myLapack.a ./

$(MY_LIBS)/myLapack.a:
	myLapackLib

$(OBJDIR)/count_trtrm.o: $(SRCDIR)/count_trtrm.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

$(OBJDIR)/count_lumm.o: $(SRCDIR)/count_lumm.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@
