chemik: main.cpp chemical.cpp stiff.cpp Makefile
	g++ -o $@ main.cpp stiff.cpp -O1 -Wall
