CC=g++
CPPFLAGS=-I.
LDFLAGS=-g
DEPS = bigraph.h utility.h abcore.h coreTree.h uf.h query.h
OBJ = bigraph.o main.o utility.o abcore.o coreTree.o uf.o query.o

%.o: %.cpp $(DEPS)
        $(CC) -std=c++1y $(LDFLAGS) -c -O3 -o $@ $< $(CPPFLAGS)

community: $(OBJ)
        $(CC) -std=c++1y $(LDFLAGS) -O3 -pthread -o $@ $^ $(CPPFLAGS)


clean :
        rm -f query *.o%