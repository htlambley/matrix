CC=gcc
CFLAGS=-lm -Wall -pg -O3

.c.o:
	$(CC) -c -o $@ $< $(CFLAGS)

all: matrix.o main.o
	$(CC) $(CFLAGS) main.o matrix.o
