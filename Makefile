CC=g++
CPPFLAGS=-I.
LDFLAGS=-g
DEPS = bigraph.h utility.h abcore.h baseline.h tabcore_baseline.h
OBJ = bigraph.o main.o utility.o abcore.o baseline.o tabcore_baseline.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y $(LDFLAGS) -c -O3 -o $@ $< $(CPPFLAGS)

tabcore: $(OBJ)
	$(CC) -std=c++1y $(LDFLAGS) -O3 -pthread -o $@ $^ $(CPPFLAGS)


clean :
	rm -f tabcore *.o