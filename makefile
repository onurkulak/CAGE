CC=g++
CFLAGS= -static -O0 -c -m64 -I/opt/ibm/ILOG/CPLEX_Studio_Community127/cplex/include/ -mcmodel=medium
CFLAGS_Link=  -O0 -m64 -I/opt/ibm/ILOG/CPLEX_Studio_Community127/cplex/include/ -L/opt/ibm/ILOG/CPLEX_Studio_Community127/cplex/lib/x86-64_linux/static_pic/ -lcplex -lm -lpthread -mcmodel=medium

cage: cage.o
	$(CC) -o cage cage.o $(CFLAGS_Link)
cage.o: cage.cpp
	$(CC) $(CFLAGS) cage.cpp
