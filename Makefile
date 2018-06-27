all: example1 example2 example3

argument_parser.o: argument_parser.cpp
	g++ -c $< -o $@

kiss_fft.o: kiss_fft.c
	gcc -c $< -o $@

reporter.o: reporter.cpp
	g++ -c $< -o $@

example1.o: example1.cpp FAI_2DFSDEsolver.h
	g++ -c $< -o $@

example1_driver.o: example1_driver.cpp
	g++ -c $< -o $@

example1: example1_driver.o argument_parser.o example1.o kiss_fft.o reporter.o
	g++ $^ -o $@

example2.o: example2.cpp FAI_2DFSDEsolver.h
	g++ -c $< -o $@

example2_driver.o: example2_driver.cpp
	g++ -c $< -o $@

example2: example2_driver.o argument_parser.o example2.o kiss_fft.o reporter.o
	g++ $^ -o $@

example3.o: example3.cpp FAI_2DFSDEsolver.h
	g++ -c $< -o $@

example3_driver.o: example3_driver.cpp
	g++ -c $< -o $@

example3: example3_driver.o argument_parser.o example3.o kiss_fft.o reporter.o
	g++ $^ -o $@

clean:
	@rm -f example1 example2 example3 *.o
