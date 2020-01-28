#############################################################################
# Makefile for 'DiskDyn'                                                    #
#                                                                           #
# Andras Gaspar                                                             #
# Steward Observatory, University of Arizona                                #
# email: agaspar@as.arizona.edu                                             #
# phone: 1-(520)-626-9898                                                   #
#                                                                           #
# Tucson, Az, US                                                            #
#############################################################################

NVCC    = nvcc

LIBS = -L ../cfitsio
INCLUDE = -I ../cfitsio

CFLAGS  = -O3 -D_GNU_SOURCE
NVFLAGS = -gencode arch=compute_35,code=sm_35 -m64 -rdc=true
LINKS   = -lcudart -lcuda -lm -lstdc++ -lcfitsio

NVCALL = $(NVCC) $(CFLAGS) $(NVFLAGS) $(LIBS) $(INCLUDE)

NVOBJECTS = ddyn.o error.o read_in.o print.o cuda.o grid.o init.o gpu_functions.o functions.o ran3.o fits.o getMemorySize.o

all: ddyn
	@if [ -d "../bin" ] ; then mv ddyn ../bin ; fi
	@echo "Compilation successful! Executable is in the bin directory."

%.o : %.cu ddyn.cuh Makefile
	$(NVCALL) -c $<

ddyn: $(NVOBJECTS) ddyn.cu ddyn.cuh
	$(NVCALL) -o ddyn $(NVOBJECTS) $(LINKS)

clean:
	rm -f *.o

nuke:
	rm -f *.o ddyn
