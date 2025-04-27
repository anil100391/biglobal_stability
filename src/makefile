CFLAGs	   = 
FFLAGS     = 
HOME       = /home/nebula/code/biglobal/src
CPPFLAGS   = -I$(HOME)/include
FPPFLAGS   =
CLEANFILES = main
LIBD       = ../lib

include $(SLEPC_DIR)/lib/slepc/conf/slepc_common

main: main.o
	${CLINKER} -o main main.o $(LIBD)/util.o -L$(LIBD)/ -lbiglobal ${SLEPC_EPS_LIB}	
	rm *.o

test: test.o
	${CLINKER} -o test test.o $(LIBD)/util.o -L$(LIBD)/ -lbiglobal ${SLEPC_EPS_LIB}	
	rm *.o

all: cheb ndr bc eig ose util
	ar q libbiglobal.a *.o
#	ar q libbiglobal.a fd_grid.o

cheb:  
	$(PCC) -c cheb.c $(CCPPFLAGS) $(PCC_FLAGS)

ndr:
	$(PCC) -c ndr.c $(CCPPFLAGS) $(PCC_FLAGS)

bc:
	$(PCC) -c bc.c $(CCPPFLAGS) $(PCC_FLAGS)
	
eig:
	$(PCC) -c eigen.c $(CCPPFLAGS) $(PCC_FLAGS)

ose:
	$(PCC) -c ose.c $(CCPPFLAGS) $(PCC_FLAGS)

util:
	$(PCC) -c memutils/util.c $(CCPPFLAGS) $(PCC_FLAGS)

install:
	mv libbiglobal.a ../lib/
	mv util.o ../lib/
	rm *.o


