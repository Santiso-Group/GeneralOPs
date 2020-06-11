
# Note: Make mymol and crystdist before making peakdist

PEAKDIST_PATH = `pwd`
# This should point to the directory containing the mymol library (mymol.a):
MYMOL_LIB_PATH = ../mymol/lib/
# This should point to the directory containing the crystdist library (crystdist.a):
CRYSTDIST_LIB_PATH = ../crystdist/lib/
# This should point to the directory containing the randomc library (randomc.a):
RANDOMC_LIB_PATH = ../randomc/lib/
# This should point to the directory containing the "common" folder:
COMMON_ROOT_PATH = "$(PEAKDIST_PATH)/../"
# This should point to the directory containing the "mymath" folder:
MYMATH_ROOT_PATH = "$(PEAKDIST_PATH)/../"
# This should point to the directory containing the "mymol" folder:
MYMOL_ROOT_PATH = "$(PEAKDIST_PATH)/../"
# This should point to the directory containing the "crystdist" folder:
CRYSTDIST_ROOT_PATH = "$(PEAKDIST_PATH)/../"

CC = g++ -m32
CFLAGS = -c -Wall -I$(COMMON_ROOT_PATH) -I$(MYMATH_ROOT_PATH) \
                  -I$(MYMOL_ROOT_PATH) -I$(CRYSTDIST_ROOT_PATH) \

all: bin/peakdist 

bin/peakdist: src/peakdist.o $(MYMOL_LIB_PATH)/mymol.a $(CRYSTDIST_LIB_PATH)/crystdist.a \
              $(RANDOMC_LIB_PATH)/randomc.a
	mkdir -p bin
	$(CC) src/peakdist.o $(MYMOL_LIB_PATH)/mymol.a $(CRYSTDIST_LIB_PATH)/crystdist.a \
        $(RANDOMC_LIB_PATH)/randomc.a -o bin/peakdist

src/peakdist.o: src/peakdist.cpp
	$(CC) $(CFLAGS) src/peakdist.cpp -o src/peakdist.o

clean:
	rm -rf src/*.o bin
