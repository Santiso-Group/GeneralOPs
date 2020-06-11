
MYMATH_PATH = `pwd`
# This should point to the directory containing the "common" folder:
COMMON_ROOT_PATH = "$(MYMATH_PATH)/../"
# This should point to the directory containing the "randomc" folder:
RANDOMC_ROOT_PATH = "$(MYMATH_PATH)/../"
MYMATH_ROOT_PATH = "$(MYMATH_PATH)/../"

CC = g++ -m32
CFLAGS = -c -Wall -I$(COMMON_ROOT_PATH) -I$(MYMATH_ROOT_PATH) \
                  -I$(RANDOMC_ROOT_PATH)

all:
	@echo "mymath is a collection of headers - no need to make."
	@echo "Use 'make debug' to generate a debug executable in the 'bin' directory."

debug: src/mymath.o
	mkdir -p bin
	$(CC) src/mymath.o $(RANDOMC_ROOT_PATH)/randomc/lib/randomc.a \
        -o bin/mymath

bootstrap: src/bootstrap.o
	mkdir -p bin
	$(CC) src/bootstrap.o $(RANDOMC_ROOT_PATH)/randomc/lib/randomc.a \
              -o bin/bootstrap

src/bootstrap.o: src/bootstrap.cpp
	$(CC) $(CFLAGS) src/bootstrap.cpp -o src/bootstrap.o

src/mymath.o: src/mymath.cpp
	$(CC) $(CFLAGS) src/mymath.cpp -o src/mymath.o

clean:
	rm -rf src/*.o bin 
	
 
