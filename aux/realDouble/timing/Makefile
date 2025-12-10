include ../make.inc
BINDIR=bin
OBJDIR=obj
SRCDIR=src

.PHONY: all clean libs tex
all: libs $(BINDIR)/time_dlarft.exe $(BINDIR)/time_dlarfb0c2.exe $(BINDIR)/time_dorgqr.exe \
	$(BINDIR)/time_dorgqr_k.exe $(BINDIR)/time_gemm.exe $(BINDIR)/time_orgkr.exe \
	#tex
#all: libs $(BINDIR)/time_dlarft.exe $(BINDIR)/time_dlarfb0c2.exe $(BINDIR)/time_dorgqr.exe $(BINDIR)/time_dorgqr_k.exe $(BINDIR)/time_gemm.exe tex

clean:
	+$(MAKE) -C referenceFiles clean
	rm $(OBJDIR)/*.o $(BINDIR)/*.exe

libs:
	+$(MAKE) -C referenceFiles all

tex:
	+$(MAKE) -C tex all

$(BINDIR)/time_dlarft.exe: $(OBJDIR)/time_dlarft.o referenceFiles/myRefLib.a myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^  $(ALL_OPT_LIB) -o $@

$(BINDIR)/time_dlarfb0c2.exe: $(OBJDIR)/time_dlarfb0c2.o referenceFiles/myRefLib.a myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^  $(ALL_OPT_LIB) -o $@

$(BINDIR)/time_dorgqr.exe: $(OBJDIR)/time_dorgqr.o referenceFiles/myRefLib.a myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^  $(ALL_OPT_LIB) -o $@

$(BINDIR)/time_dorgqr_k.exe: $(OBJDIR)/time_dorgqr_k.o referenceFiles/myRefLib.a myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^  $(ALL_OPT_LIB) -o $@

$(BINDIR)/time_orgkr.exe: $(OBJDIR)/time_orgkr.o referenceFiles/myRefLib.a myLapack.a
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^  $(ALL_OPT_LIB) -o $@

$(BINDIR)/time_gemm.exe: $(OBJDIR)/time_gemm.o
	$(CPP) $(CFLAGS) $(ALL_OPT_INC) -lgfortran $^ $(ALL_OPT_LIB) -o $@

myLapack.a: $(MY_LIBS)/myLapack.a
	cp $(MY_LIBS)/myLapack.a ./
	ar d myLapack.a dlarft.o zlarft.o clarft.o slarft.o dorgqr.o sorgqr.o dorglq.o sorglq.o

$(MY_LIBS)/myLapack.a:
	myLapackLib

referenceFiles/myRefLib.a:
	+$(MAKE) -C referenceFiles all

$(OBJDIR)/time_dlarft.o: $(SRCDIR)/time_dlarft.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

$(OBJDIR)/time_dlarfb0c2.o: $(SRCDIR)/time_dlarfb0c2.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

$(OBJDIR)/time_dorgqr.o: $(SRCDIR)/time_dorgqr.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

$(OBJDIR)/time_dorgqr_k.o: $(SRCDIR)/time_dorgqr_k.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

$(OBJDIR)/time_orgkr.o: $(SRCDIR)/time_orgkr.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@

$(OBJDIR)/time_gemm.o: $(SRCDIR)/time_gemm.f
	$(FC) -c $(FFLAGS) -O0 $< -o $@
