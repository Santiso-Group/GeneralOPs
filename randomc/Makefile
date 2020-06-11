
RANDOMC_PATH = `pwd`
RANDOMC_ROOT_PATH = "$(RANDOMC_PATH)/../"

CC = g++ -m32
CFLAGS = -c -Wall -I$(RANDOMC_ROOT_PATH)

all: lib/randomc.a

lib/randomc.a: src/mother.o src/mersenne.o src/rancombi.o
	mkdir -p lib
	ar cq lib/randomc.a src/mother.o src/mersenne.o src/rancombi.o

src/mother.o: src/mother.cpp
	$(CC) $(CFLAGS) src/mother.cpp -o src/mother.o

src/mersenne.o: src/mersenne.cpp
	$(CC) $(CFLAGS) src/mersenne.cpp -o src/mersenne.o

src/rancombi.o: src/rancombi.cpp
	$(CC) $(CFLAGS) src/rancombi.cpp -o src/rancombi.o

clean:
	rm -rf src/*o $(DATA_INCLUDE_FILE) lib bin
	
 
