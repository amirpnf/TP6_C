# Makefile for compiling collision.c into a shared library

CC =gcc
CFLAGS =-c -fPIC -Wall -ansi -O2
LDFLAGS =-shared
TARGET =libcollision.so

all: $(TARGET)

$(TARGET): collision.o
	$(CC) $(LDFLAGS) -o $@ $^

collision.o: collision.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f *.o $(TARGET)
