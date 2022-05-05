CC=g++
CPPFLAGS=-I.
LDFLAGS=-g
DEPS = bigraph/bigraph.h utility/utility.h abcore/abcore.h baseline/baseline.h baseline/tabcore_baseline.h
OBJ = bigraph/bigraph.o main.o utility/utility.o abcore/abcore.o baseline/baseline.o baseline/tabcore_baseline.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y $(LDFLAGS) -c -O3 -o $@ $< $(CPPFLAGS)

tabcore: $(OBJ)
	$(CC) -std=c++1y $(LDFLAGS) -O3 -pthread -o $@ $^ $(CPPFLAGS)

clean :
	rm main.o
	rm bigraph/*.o
	rm utility/*.o
	rm abcore/*.o
	rm baseline/*.o
	rm tabcore