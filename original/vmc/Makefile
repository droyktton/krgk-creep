OBJS = system.o mainalbe.o mt.o minimize.o vmc.o actions.o

G2DIR=/usr/

SRCS=${OBJS:.o=.cpp}
 
CFLAGS= -I$(G2DIR) -I$(BZDIR) -I/usr/local/include -I/usr/X11R6/include #-pg -g

DEFS=  -DOPTIMIZE2 -DOPTIMIZE3 -DEXACT -DOPTIMIZE5 -DINSTANTONIC_WALK #-DFIXEDTM
#-DSTRONG_PINNING
#-DPERIODIC_CHAIN -DPERIODIC_DISORDER # #-DPINNING_CENTER  -DRANDOM_FIELD 

DEFS2= #-DLOCAL_MOVES -DINSTANTONIC_WALK  

CFLAGS2= -axW -xW -ipo -ansi -Xc -O3 #-ftemplate-depth-30 
 
LIBS=  -lblitz -lm -lfftw3 -O3  

CC= g++

# to see the movie uncomment 3 following lines
#CC= g++-4.4 #-DGRAPHICS
GLIBS= -lg2 -lX11 -lgd #-L/usr/lib/x86_64-linux-gnu/
GDEF = -DGRAPHICS #-DMOVIE_PNG #-DEPS #-DCORRELATED_DISORDER

all:   $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) $(GLIBS) -o VMC

system.o : vmc.hpp mt.cpp system.cpp 
	$(CC) -c $(DEFS)  $(GDEF) $(CFLAGS)  system.cpp

minimize.o : vmc.hpp minimize.cpp
	$(CC) -c $(DEFS)  $(GDEF) $(CFLAGS)  minimize.cpp

vmc.o : vmc.hpp vmc.cpp
	$(CC) -c $(DEFS)  $(GDEF) $(CFLAGS)  vmc.cpp

mainalbe.o : vmc.hpp mainalbe.cpp
	$(CC) -c $(DEFS) $(GDEF) $(CFLAGS)  mainalbe.cpp

actions.o : vmc.hpp actions.cpp
	$(CC) -c $(DEFS) $(GDEF) $(CFLAGS)  actions.cpp

mt.o : mt.cpp 
	$(CC) -c $(DEFS) $(CFLAGS)  mt.cpp

clean:
	rm -f $(OBJS) *.o VMC

print:
	a2ps *.hpp $(SRCS) Makefile
        
archive:
	tar -cf vmc.tar *pp 

edata:
	rm *.dat zzz *.fft
