
# Note: Make mymol and crystdist before making getpeaks

GETPEAKS_PATH = `pwd`
# This should point to the directory containing the mymol library (mymol.a):
MYMOL_LIB_PATH = ../mymol/lib/
# This should point to the directory containing the crystdist library (crystdist.a):
CRYSTDIST_LIB_PATH = ../crystdist/lib/
# This should point to the directory containing the "common" folder:
COMMON_ROOT_PATH = "$(GETPEAKS_PATH)/../"
# This should point to the directory containing the "mymath" folder:
MYMATH_ROOT_PATH = "$(GETPEAKS_PATH)/../"
# This should point to the directory containing the "mymol" folder:
MYMOL_ROOT_PATH = "$(GETPEAKS_PATH)/../"
# This should point to the directory containing the "crystdist" folder:
CRYSTDIST_ROOT_PATH = "$(GETPEAKS_PATH)/../"

CC = g++ -m32
CFLAGS = -c -Wall -I$(COMMON_ROOT_PATH) -I$(MYMATH_ROOT_PATH) \
                  -I$(MYMOL_ROOT_PATH) -I$(CRYSTDIST_ROOT_PATH) \
                  -I$(GETPEAKS_PATH)

all: bin/getpeaks 

bin/getpeaks: src/getpeaks.o $(MYMOL_LIB_PATH)/mymol.a $(CRYSTDIST_LIB_PATH)/crystdist.a
	mkdir -p bin
	$(CC) src/getpeaks.o $(MYMOL_LIB_PATH)/mymol.a $(CRYSTDIST_LIB_PATH)/crystdist.a -o bin/getpeaks

src/getpeaks.o: src/getpeaks.cpp
	$(CC) $(CFLAGS) src/getpeaks.cpp -o src/getpeaks.o

clean:
	rm -rf src/*.o bin
