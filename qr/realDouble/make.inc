# This is a path to where the LAPACK release library exists. Below is mine
# I set an environment variable that I use. You can do the same or type out the
# whole path

LAPACKROOT=$(CUDENVERPATH)/LAPACK/lapack-3.11.0

LAPACKEINC=$(LAPACKROOT)/LAPACKE/include
CBLASINC=$(LAPACKROOT)/CBLAS/include

LAPACKLIB=$(LAPACKROOT)/liblapack.a
LAPACKELIB=$(LAPACKROOT)/liblapacke.a
BLASLIB=$(LAPACKROOT)/librefblas.a
CBLASLIB=$(LAPACKROOT)/libcblas.a

CC=gcc
FC=gfortran

CFLAGS=-g
FFLAGS=-fimplicit-none -g

# This is the compilation flags needed to compile with optimized blas
# For intel machines, this whatever is needed for MKL.
#	  On Linux, this is done by instally mkl (typically through the package manager) and then running 
#	  ```mkl_link_tool gfortran $(CFLAGS) my_dorgqr.o test.o my_dlarfb.o```
#	  and then taking the ending bit after 'my_dlarfb.o'
# For AMD machines, this is whatever is needed for AOCL
#         You can be done something like
#         OPTBLAS=<path-to-libFLAME-library>/libflame.so <path-to-BLIS-library>/libblis.so
#         or to follow the up to date recommendation from the AOCL user guide in the chapter discussing using 
#         AOCL-libFLAME
#OPTBLAS=-m64  -I"/opt/intel/oneapi/mkl/2023.2.0/include"  -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#OPTBLAS=$(AOCLROOT)/lib/libflame.so $(AOCLROOT)/lib/libblis.so
OPTBLAS=/home/lidesia/Documents/CUDenver/RAShip/gitRepos/aocl-utils/build/libaoclutils.a /home/lidesia/Documents/CUDenver/RAShip/gitRepos/libflame/lib/x86_64-unknown-linux-gnu/libflame.a /home/lidesia/Documents/CUDenver/RAShip/gitRepos/blis/lib/amd64/libblis.a
