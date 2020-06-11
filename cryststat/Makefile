
# Note: Make mymol and crystdist before making cryststat

CRYSTSTAT_PATH = `pwd`
# This should point to the directory containing the mymol library (mymol.a):
MYMOL_LIB_PATH = ../mymol/lib/
# This should point to the directory containing the crystdist library (crystdist.a):
CRYSTDIST_LIB_PATH = ../crystdist/lib/
# This should point to the directory containing the "common" folder:
COMMON_ROOT_PATH = "$(CRYSTSTAT_PATH)/../"
# This should point to the directory containing the "mymath" folder:
MYMATH_ROOT_PATH = "$(CRYSTSTAT_PATH)/../"
# This should point to the directory containing the "mymol" folder:
MYMOL_ROOT_PATH = "$(CRYSTSTAT_PATH)/../"
# This should point to the directory containing the "crystdist" folder:
CRYSTDIST_ROOT_PATH = "$(CRYSTSTAT_PATH)/../"

CC = g++ -m32
CFLAGS = -c -Wall -I$(COMMON_ROOT_PATH) -I$(MYMATH_ROOT_PATH) \
                  -I$(MYMOL_ROOT_PATH) -I$(CRYSTDIST_ROOT_PATH) 

all: bin/cryststat 

bin/cryststat: src/cryststat.o $(MYMOL_LIB_PATH)/mymol.a $(CRYSTDIST_LIB_PATH)/crystdist.a
	mkdir -p bin
	$(CC) src/cryststat.o $(MYMOL_LIB_PATH)/mymol.a $(CRYSTDIST_LIB_PATH)/crystdist.a -o bin/cryststat

src/cryststat.o: src/cryststat.cpp
	$(CC) $(CFLAGS) src/cryststat.cpp -o src/cryststat.o

clean:
	rm -rf src/*.o bin
