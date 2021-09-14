all: mim

mim: *.cpp imm/*.cpp imm/*.h
	g++ -O3 -std=c++0x -DDISCRETE mim.cpp imm/sfmt/SFMT.c -o mim -lpthread
test: *.cpp imm/*.cpp imm/*.h
	g++ -O0 -g -std=c++0x -DDISCRETE mim.cpp imm/sfmt/SFMT.c -o mim -lpthread
synth: *.cpp
	g++ -O3 -std=c++11 synth_multi.cpp -o synth -lpthread
