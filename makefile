
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic


cpps=$(shell ls source/*.cpp) 
headers=$(shell ls include/*.h)
sources=$(cpps) $(headers)
includes=-I.
 
CCC = g++
CC  = gcc  
CCOPT = -fPIC -fexceptions -DNDEBUG -DIL_STD -w -g -Wall -Wextra -ldl
COPT  = -fPIC

CPLEXDIR      = /home/joe/CPLEX_Studio126/cplex
CONCERTDIR    = /home/joe/CPLEX_Studio126/concert

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

LEMONLIBDIR   = /usr/local/lib
LEMONINCDIR   = /usr/local/include 

BOOST = /home/joe/boost_1_78_0

LEMONFLAGS	= -I$(LEMONINCDIR) -L$(LEMONLIBDIR) -I$(BOOST) -lemon

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -m64 -lm -lpthread -ldl

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(LEMONINCDIR) $(LEMONFLAGS)
	


all:
#/*ajouter les fichiers (FLAGS)*/ 
	$(CCC) -c -g $(CCFLAGS) $(cpps) $(includes)
	$(CCC) $(CCFLAGS) *.o -g -o Survivability $(CCLNFLAGS)
	rm -rf *.o *~ ^
#-lemon

clean:
	rm ./obj/*.o ag *.oo core* *.o *.sav *.lp -f
