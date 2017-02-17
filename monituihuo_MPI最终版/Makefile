#!/bin/bash

CC=mpiicc  
CFLAGS=  -O2  -g 
LDFLAGS=  -O2  -g
#CFLAGS=     
#LDFLAGS=  
#CFLAGS= -openmp  -guide -vec-report=3 -ansi -O3 -Wall
#LDFLAGS= -openmp  -guide  -vec-report=3 -ansi -lm -Wall
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
