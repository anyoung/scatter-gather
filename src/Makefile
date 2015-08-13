CC=gcc
CFLAGS=-g -fPIC -Wall

OBJS=scatgat.o sg_access.o

.PHONY: all clean

all: libscatgat.so

clean:
	rm -f *.o
	rm -f *.so

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

libscatgat.so: $(OBJS)
	$(CC) -shared -o $@ $^
