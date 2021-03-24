
MYMOL_PATH = `pwd`
# This should point to the directory containing the "common" folder:
COMMON_ROOT_PATH = "$(MYMOL_PATH)/../"
# This should point to the directory containing the "mymath" folder:
MYMATH_ROOT_PATH = "$(MYMOL_PATH)/../"
MYMOL_ROOT_PATH = "$(MYMOL_PATH)/../"
DATA_INCLUDE_FILE = "$(MYMOL_PATH)/include/datafiles.h"

# Compile in 32-bit if having trouble reading dcd files
#CC = g++ -m32 -fast
#CC = g++ -fast
#CC = clang++ -m32
CC = g++ -m32
#CC = clang++ -O4 -std=c++11 -stdlib=libc++
CFLAGS = -c -Wall -I$(COMMON_ROOT_PATH) -I$(MYMATH_ROOT_PATH) -I$(MYMOL_ROOT_PATH)

all: lib/mymol.a

lib/mymol.a: src/atom.o src/dcdfile.o src/geometry.o src/lattice.o \
             src/mol2file.o src/molecule.o src/pdbfile.o src/psffile.o \
             src/configfile.o src/pairpotential.o src/fieldfile.o src/historyfile.o \
             src/discretefield.o src/grid.o src/vtkfile.o src/xyzfile.o
	mkdir -p lib
	ar cq lib/mymol.a src/atom.o src/dcdfile.o src/geometry.o src/lattice.o \
                          src/mol2file.o src/molecule.o src/pdbfile.o src/psffile.o \
                          src/pairpotential.o src/configfile.o src/fieldfile.o src/historyfile.o \
                          src/discretefield.o src/grid.o src/vtkfile.o src/xyzfile.o

debug: lib/mymol.a src/mymol.o 
	mkdir -p bin
	$(CC) src/mymol.o lib/mymol.a -o bin/mymol

src/mymol.o: src/mymol.cpp
	$(CC) $(CFLAGS) src/mymol.cpp -o src/mymol.o

src/atom.o: src/atom.cpp 
	rm -f $(DATA_INCLUDE_FILE)
	echo "std::string const ATOM_DATABASE = \"$(MYMOL_PATH)/dat/atom.dat\";" >> $(DATA_INCLUDE_FILE)
	echo >> $(DATA_INCLUDE_FILE)
	$(CC) $(CFLAGS) src/atom.cpp -o src/atom.o

src/dcdfile.o: src/dcdfile.cpp
	$(CC) $(CFLAGS) src/dcdfile.cpp -o src/dcdfile.o

src/geometry.o: src/geometry.cpp
	$(CC) $(CFLAGS) src/geometry.cpp -o src/geometry.o

src/lattice.o: src/lattice.cpp
	$(CC) $(CFLAGS) src/lattice.cpp -o src/lattice.o

src/mol2file.o: src/mol2file.cpp
	$(CC) $(CFLAGS) src/mol2file.cpp -o src/mol2file.o

src/molecule.o: src/molecule.cpp
	$(CC) $(CFLAGS) src/molecule.cpp -o src/molecule.o

src/pairpotential.o: src/pairpotential.cpp
	$(CC) $(CFLAGS) src/pairpotential.cpp -o src/pairpotential.o

src/pdbfile.o: src/pdbfile.cpp
	$(CC) $(CFLAGS) src/pdbfile.cpp -o src/pdbfile.o

src/psffile.o: src/psffile.cpp
	$(CC) $(CFLAGS) src/psffile.cpp -o src/psffile.o

src/configfile.o: src/configfile.cpp
	$(CC) $(CFLAGS) src/configfile.cpp -o src/configfile.o

src/fieldfile.o: src/fieldfile.cpp
	$(CC) $(CFLAGS) src/fieldfile.cpp -o src/fieldfile.o

src/historyfile.o: src/historyfile.cpp
	$(CC) $(CFLAGS) src/historyfile.cpp -o src/historyfile.o

src/discretefield.o: src/discretefield.cpp
	$(CC) $(CFLAGS) src/discretefield.cpp -o src/discretefield.o

src/grid.o: src/grid.cpp
	$(CC) $(CFLAGS) src/grid.cpp -o src/grid.o

src/vtkfile.o: src/vtkfile.cpp
	$(CC) $(CFLAGS) src/vtkfile.cpp -o src/vtkfile.o

src/xyzfile.o: src/xyzfile.cpp
	$(CC) $(CFLAGS) src/xyzfile.cpp -o src/xyzfile.o

clean:
	rm -rf src/*o $(DATA_INCLUDE_FILE) lib bin
	
 
