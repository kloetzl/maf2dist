CPPFLAGS += -Wall -Wextra -std=gnu++14
CXXFLAGS += -O3 -ggdb -march=native


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

