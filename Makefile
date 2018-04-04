all: example2 example3

example2: example2.cpp
	g++ example2.cpp kiss_fft.c -o example2

example3: example3.cpp
	g++ example3.cpp kiss_fft.c -o example3
