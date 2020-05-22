CC = gcc
MCC = /usr/lib64/openmpi/bin/mpicc
CFLAGS = -l /usr/lib64/libgsl.so /usr/lib64/libgslcblas.so
INCLUDE_DIR = -I /home/zongboxu/MCMS/Density
GFLAGS = -Wall -O2 
OBJS = /home/zongboxu/MCMS/Density/procdfcs.c
#-I/home/lmargerin/lib/include 

mcpray : mcpray.c
	$(MCC) $(GFLAGS) -o mcpray ./Lambda/mcpray.c -lm -lgslcblas -lgsl
mcprayfs : mcprayfs.c
	$(MCC) $(GFLAGS) -o mcprayfs ./Lambda/mcprayfs.c -lm -lgslcblas -lgsl
mcprayden : ./Density/mcprayavden.c
	$(MCC) $(GFLAGS) ${INCLUDE_DIR} -o mcprayden ./Density/mcprayavden.c ${OBJS} -lm -lgslcblas -lgsl
mcbrayden : ./Density/mcbrayavden.c
	$(MCC) $(GFLAGS) ${INCLUDE_DIR} -o mcbrayden ./Density/mcbrayavden.c ${OBJS} -lm -lgslcblas -lgsl
mcbraydenfs : ./Density/mcbrayavdenfs.c
	$(MCC) $(GFLAGS) ${INCLUDE_DIR} -o mcbraydenfs ./Density/mcbrayavdenfs.c ${OBJS} -lm -lgslcblas -lgsl
mcbraydenlb : ./Density/mcbrayavdenlb.c
	$(MCC) $(GFLAGS) ${INCLUDE_DIR} -o mcbraydenlb ./Density/mcbrayavdenlb.c ${OBJS} -lm -lgslcblas -lgsl