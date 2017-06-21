CPPFLAGS += -Wall -Wextra -std=gnu++14
CXXFLAGS += -O2 -ggdb -march=native -ftree-vectorize


.PHONY: all clean format

all: maf2dist

%.o: %.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $^

maf2dist: maf2dist.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f maf2dist *.o

format:
	clang-format -i *.cxx

