CC=clang
CFLAGS=-lm -Wall -Wpedantic -pg -O3 -ferror-limit=1

.c.o:
	$(CC) -c -o $@ $< $(CFLAGS)

all: matrix.o main.o
	$(CC) $(CFLAGS) main.o matrix.o
