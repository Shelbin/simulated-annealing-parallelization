#!/bin/bash

CC=icc  
CFLAGS=  -O2  -g 
LDFLAGS=  -O2  -g

EXEC=monituihuo
OBJ1=monituihuo.o 

all: $(EXEC)

monituihuo : $(OBJ1)
	$(CC) -o $@ $^ $(LDFLAGS)

 
##########################################
# Generic rules
##########################################

%.o: %.c %.h
	$(CC) -o $@ -c $< $(CFLAGS)
	
%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)
	
clean:
	rm -f *.o *~ $(EXEC)
