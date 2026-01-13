include ../make.inc
OBJDIR=obj
SRCDIR=src
BINDIR=bin

SRCFILES=$(wildcard $(SRCDIR)/*.f)
OBJFILES=$(patsubst $(SRCDIR)/%.f,$(OBJDIR)/%.o,$(SRCFILES))
EXEFILES=$(patsubst $(SRCDIR)/%.f,$(BINDIR)/%.exe,$(SRCFILES))

.PHONY: all clean
#all: $(BINDIR)/test_dstrk.exe $(BINDIR)/test_dst3rk.exe $(BINDIR)/test_dtrmvoop.exe $(BINDIR)/test_dtrmmoop_lvl2.exe $(BINDIR)/test_dorgkr.exe
all: $(EXEFILES)

clean:
	rm $(OBJDIR)/*.o $(BINDIR)/*.exe

cleanlib:
	rm myLapack.a

myLapack.a: $(MY_LIBS)/myLapack.a
	cp $(MY_LIBS)/myLapack.a ./

$(BINDIR)/%.exe: $(OBJDIR)/%.o myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(ALL_OPT_LIB) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

