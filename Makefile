CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
BINDIR=/usr/local/bin
objects = main.o liheng.o statfq.o NoNregion.o
f-tools:$(objects)
	$(CC) $(CFLAGS) -o f-tools $(objects) -lz -lm
	rm -rf *.o
statfq.o:statfq.c kseq.h
	$(CC) $(CFLAGS) -c statfq.c -lz -lm
NoNregion.o:NoNregion.c
	$(CC) $(CFLAGS) -c NoNregion.c -lm
liheng.o:liheng.c khash.h kseq.h
	$(CC) $(CFLAGS) -c liheng.c -o $@ -lz -lm
main.o:main.c main.h
	$(CC) $(CFLAGS) -c main.c -lz -lm
clean:
	rm -fr gmon.out *.o ext/*.o trimadap *~ *.a *.dSYM session*
