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

CFLAGS  = -std=c99 -O3 -D_GNU_SOURCE
LINKS   = -lm

all: runmie

runmie: runmie.c mie_single.c Makefile
	$(CC) $(CFLAGS) runmie.c mie_single.c -o runmie $(LINKS)

clean:
	rm -f *.o

nuke:
	rm -f *.o runmie
