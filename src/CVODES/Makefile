# Makefile to compile libcvodes.a
# author: TC
#
#
CC = gcc
CFLAGS = -g -O2 -fPIC
OUT = libcvodes.a


SRC := $(wildcard *.c)
#SRC = cvodea.c	cvodea_io.c 	cvodes_band.c cvodes_bandpre.c 	cvodes_bbdpre.c cvodes_dense.c 	cvodes_diag.c\
cvodes_io.c	cvodes_spbcgs.c	cvodes_spgmr.c	cvodes_spils.c	cvodes_sptfqmr.c	fnvector_serial.c\
	ida_band.c	ida_bbdpre.c	ida_dense.c	ida_ic.c	ida_io.c	ida_spbcgs.c	ida_spgmr.c	ida_spils.c	ida_sptfqmr.c\
	nvector_serial.c	sundials_band.c	sundials_cvodes.c	sundials_dense.c	sundials_ida.c	sundials_iterative.c	sundials_math.c\
	sundials_nvector.c	sundials_smalldense.c	sundials_spbcgs.c	sundials_spgmr.c	sundials_sptfqmr.c

OBJ = $(SRC:.c=.o)

.SUFFIXES : .c



.c.o:
	$(CC) $(CFLAGS) -c $< -o $@


all: $(OUT) sharelib

$(OUT): $(OBJ)
	ar rcs $(OUT) $(OBJ)

sharelib:
	$(CC) -fPIC -shared $(OBJ) -o libcvodes.so

	
clean:
	rm -f $(OBJ) $(OUT)

