CFLAGs	   = 
FFLAGS     = 
HOME       = /home/nebula/Dropbox/ME_project/biglobal/src
CPPFLAGS   = -I$(HOME)/include
FPPFLAGS   =
CLEANFILES = main
LIBD       = /home/nebula/Dropbox/ME_project/biglobal/lib

include $(SLEPC_DIR)/conf/slepc_common

main_ose: main_ose.o chkopts
	${CLINKER} -o main_ose main_ose.o $(LIBD)/util.o -L$(LIBD)/ -lbiglobal ${SLEPC_EPS_LIB}	
	rm *.o

try: try.o chkopts
	${CLINKER} -o try try.o $(LIBD)/util.o -L$(LIBD)/ -lbiglobal ${SLEPC_EPS_LIB}	
	rm *.o

coflow_wakes: coflow_wakes.o chkopts
	${CLINKER} -o coflow_wakes coflow_wakes.o $(LIBD)/util.o -L$(LIBD)/ -lbiglobal ${SLEPC_EPS_LIB}
	rm *.o

check_diff: check_diff.o chkopts
	${CLINKER} -o check_diff check_diff.o $(LIBD)/util.o -L$(LIBD)/ -lbiglobal ${SLEPC_EPS_LIB}
	rm *.o
