CC=gcc
CXX=g++
CFLAGS=-O3 -Wall -fPIC -g
CLIBS=-lm -lcfitsio -lfftw3_threads -lfftw3
LAPACK=-lopenblas -lpthread -lgfortran

WCSI=-I/usr/include/wcslib
WCSL=-L/usr/lib
WCSLIBS=-lwcs

# glib
GLIBI=-I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include -I/usr/lib/x86_64-linux-gnu/glib-2.0/include/

INCLUDES=-I. -I/usr/include $(WCSI) $(GLIBI)
LIBPATH=-L/usr/lib $(WCSL) $(LMLIB) $(L1LIB)

RANLIB=ranlib
OBJECTS=mmio.o fftutil.o myblasf.o pswf.o wkernel2.o lmwcs.o  clmfit_nocuda.o  myblas2.o  imagetf.o ptree.o flagger.o  
default:liblargefft.a
mmio.o:mmio.c largefft.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
fftutil.o:fftutil.c largefft.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
myblasf.o:myblasf.c largefft.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
myblas2.o:myblas2.c largefft.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
pswf.o:pswf.c largefft.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
wkernel2.o:wkernel2.c largefft.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
lmwcs.o:lmwcs.c largefft.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
imagetf.o:imagetf.c largefft.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
clmfit_nocuda.o:clmfit_nocuda.c largefft.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
ptree.o:ptree.c gridder.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
flagger.o:flagger.c gridder.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

liblargefft.a:$(OBJECTS)
	ar rv liblargefft.a $(OBJECTS); \
	$(RANLIB) liblargefft.a;

clean:
	rm -f *.a *.o
