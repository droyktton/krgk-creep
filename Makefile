CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -isystem third_party

build/krgk_creep: src/krgk_creep.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -o $@ $<

test: build/krgk_creep
	tests/regression.sh

clean:
	rm -rf build

.PHONY: test clean
