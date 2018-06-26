all: example2 example3

example2: example2.cpp FAI_2DFSDEsolver.h
	g++ example2.cpp kiss_fft.c -o example2

example3: example3.cpp FAI_2DFSDEsolver.h
	g++ example3.cpp kiss_fft.c -o example3

clean:
	@rm -f example2 example3
