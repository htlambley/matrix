CC=clang
LDFLAGS=-lm
CFLAGS=-Wall -Wpedantic -O3 -ferror-limit=1
SRC = $(wildcard *.c)
OBJECTS = $(SRC:.c=.o)

%o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

main: $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) $(OBJECTS) -o $@
all: main;

clean:
	rm -rf $(OBJECTS) ./a.out
