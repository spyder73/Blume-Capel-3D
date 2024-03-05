COMPILER=g++ -c
COMPILER_FLAGS=-std=c++14 -Wall -Wextra -pedantic -O3
LINKER=g++

TARGETS=main.o statistics.o blume.o triple.o random.o

all: a.out

a.out: $(TARGETS)
	$(LINKER) $^ -o $@

%.o: %.cpp
	$(COMPILER) $(COMPILER_FLAGS) -c $< -o $@

clean:
	rm -f *.o a.out

.PHONY: all clean