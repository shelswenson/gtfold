INCLUDE_PATH =
LIBRARY_PATH =
CFLAGS = $(INCLUDE_PATH) -g -Wall
LDFLAGS =  -g -Wall
CC = gcc
CPP = g++
OBJS = loader.o	LoopScoring.o RNAScoring.o TreeScoring.o StructureReader.o
LINKOBJ  = loader.o LoopScoring.o StructureReader.o TreeScoring.o
SRC = RNAscoring.cc StructureReader.c TreeScoring.c LoopScoring.c loader.cc

LIBS =

#all: RNAscoring
#RNAScoring.o: RNAScoring.cc
#	$(CPP) -c RNAScoring.cc 

#%.o: %.c
#	$(CC) -c $(CFLAGS) $< -o $@

#%.o: %.cc
#	$(CPP) -c $< -o $@

#RNAscoring: $(OBJS)
#	$(CPP)  -o $@ $^ $(LIBRARY_PATH) $(LIBS)

RNAScoring: $(LINKOBJ) RNAScoring.cc
	$(CPP) RNAScoring.cc $(LINKOBJ) -o RNAScoring


loader.o: loader.cc
	$(CPP) -c loader.cc -o loader.o $

LoopScoring.o: LoopScoring.c
	$(CC) -c LoopScoring.c -o LoopScoring.o $(CFLAGS)

StructureReader.o: StructureReader.c
	$(CC) -c StructureReader.c -o StructureReader.o $(CFLAGS)

TreeScoring.o: TreeScoring.c
	$(CC) -c TreeScoring.c -o TreeScoring.o $(CFLAGS)

clean:
	rm -f $(OBJS)



#depends:
#	makedepend $(INCLUDE_PATH) $(SRC)

# DO NOT DELETE

