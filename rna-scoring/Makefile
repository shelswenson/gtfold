INCLUDE_PATH =
LIBRARY_PATH =
CFLAGS = $(INCLUDE_PATH) -arch i386 -g -Wall
LDFLAGS = -arch i386 -g -Wall
CC = gcc

SRC = RNAScoring.c

LIBS =

OBJS = $(SRC:.cpp=.o)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

RNAScoring: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBRARY_PATH) $(LIBS)

clean:
	rm -f $(OBJS)

depends:
	makedepend $(INCLUDE_PATH) $(SRC)

# DO NOT DELETE

