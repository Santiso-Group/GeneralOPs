
# Note: Make mymol before making crystdist

CRYSTDIST_PATH = `pwd`
# This should point to the directory containing the "common" folder:
COMMON_ROOT_PATH = "$(CRYSTDIST_PATH)/../"
# This should point to the directory containing the "mymath" folder:
MYMATH_ROOT_PATH = "$(CRYSTDIST_PATH)/../"
# This should point to the directory containing the "mymol" folder: 
MYMOL_ROOT_PATH = "$(CRYSTDIST_PATH)/../"
# This should point to the directory containing the "randomc" folder:
RANDOMC_ROOT_PATH = "$(CRYSTDIST_PATH)/../"
CRYSTDIST_ROOT_PATH = "$(CRYSTDIST_PATH)/../"

CC = g++ -m32
CFLAGS = -c -Wall -I$(COMMON_ROOT_PATH) -I$(MYMATH_ROOT_PATH) -I$(MYMOL_ROOT_PATH) -I$(CRYSTDIST_ROOT_PATH) -I$(RANDOMC_ROOT_PATH)

all: lib/crystdist.a

lib/crystdist.a: src/copfile.o src/globalhistogram.o src/mmpfile.o src/opcalculator.o \
                 src/peakaccumulator.o src/peakstatistics.o src/peakclustering.o src/peaktable.o \
                 src/pktfile.o src/pointmolecule.o src/ptmfile.o src/relativeconfiguration.o \
                 src/xtpfile.o 
	mkdir -p lib
	ar cq lib/crystdist.a src/copfile.o src/globalhistogram.o src/mmpfile.o \
                              src/opcalculator.o src/peakaccumulator.o \
                              src/peakstatistics.o src/peakclustering.o \
                              src/peaktable.o src/pktfile.o src/pointmolecule.o \
                              src/ptmfile.o src/relativeconfiguration.o src/xtpfile.o

debug: lib/crystdist.a src/crystdist.o 
	mkdir -p bin
	$(CC) src/crystdist.o lib/crystdist.a $(MYMOL_ROOT_PATH)/mymol/lib/mymol.a \
	      $(RANDOMC_ROOT_PATH)/randomc/lib/randomc.a -o bin/crystdist

src/crystdist.o: src/crystdist.cpp
	$(CC) $(CFLAGS) src/crystdist.cpp -o src/crystdist.o

src/copfile.o: src/copfile.cpp
	$(CC) $(CFLAGS) src/copfile.cpp -o src/copfile.o

src/globalhistogram.o: src/globalhistogram.cpp
	$(CC) $(CFLAGS) src/globalhistogram.cpp -o src/globalhistogram.o

src/mmpfile.o: src/mmpfile.cpp
	$(CC) $(CFLAGS) src/mmpfile.cpp -o src/mmpfile.o

src/opcalculator.o: src/opcalculator.cpp
	$(CC) $(CFLAGS) src/opcalculator.cpp -o src/opcalculator.o

src/peakaccumulator.o: src/peakaccumulator.cpp
	$(CC) $(CFLAGS) src/peakaccumulator.cpp -o src/peakaccumulator.o

src/peakstatistics.o: src/peakstatistics.cpp include/peakstatistics.h
	$(CC) $(CFLAGS) src/peakstatistics.cpp -o src/peakstatistics.o

src/peakclustering.o: src/peakclustering.cpp include/peakclustering.h
	$(CC) $(CFLAGS) src/peakclustering.cpp -o src/peakclustering.o

src/peaktable.o: src/peaktable.cpp
	$(CC) $(CFLAGS) src/peaktable.cpp -o src/peaktable.o

src/pktfile.o: src/pktfile.cpp
	$(CC) $(CFLAGS) src/pktfile.cpp -o src/pktfile.o

src/pointmolecule.o: src/pointmolecule.cpp
	$(CC) $(CFLAGS) src/pointmolecule.cpp -o src/pointmolecule.o

src/ptmfile.o: src/ptmfile.cpp
	$(CC) $(CFLAGS) src/ptmfile.cpp -o src/ptmfile.o

src/relativeconfiguration.o: src/relativeconfiguration.cpp
	$(CC) $(CFLAGS) src/relativeconfiguration.cpp -o src/relativeconfiguration.o

src/xtpfile.o: src/xtpfile.cpp
	$(CC) $(CFLAGS) src/xtpfile.cpp -o src/xtpfile.o

clean:
	rm -rf *.x src/*o lib bin

