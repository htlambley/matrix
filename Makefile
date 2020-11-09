CC=clang
LDFLAGS=-lm
CFLAGS=-Wall -Wpedantic -O3 -ferror-limit=1
SRC = $(wildcard *.c)
OBJECTS = $(SRC:.c=.o)

%o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: test/main.o package
	$(CC) $(LDFLAGS) $(CFLAGS) test/main.o libmatrix.a -o test/main

package: matrix.o
	ar -cq libmatrix.a matrix.o

clean:
	rm -rf $(OBJECTS) *.a ./a.out ./test/*.o
