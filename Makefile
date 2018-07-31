all: example1 example1_BDADI example2 example2_BFSMGM example3

CC := g++
C := gcc

CXXFLAGS := -g -O3 --std=c++11

argument_parser.o: argument_parser.cpp
	${CC} ${CXXFLAGS} -c $< -o $@

kiss_fft.o: kiss_fft.c
	${C} ${CXXFLAGS} -c $< -o $@

reporter.o: reporter.cpp
	${CC} ${CXXFLAGS} -c $< -o $@

example1.o: example1.cpp FAI_2DFSDEsolver.h
	${CC} ${CXXFLAGS} -c $< -o $@

example1_driver.o: example1_driver.cpp
	${CC} ${CXXFLAGS} -c $< -o $@

example1: example1_driver.o argument_parser.o example1.o kiss_fft.o reporter.o
	${CC} ${CXXFLAGS} $^ -o $@

example1_BDADI.o: example1_BDADI.cpp
	${CC} ${CXXFLAGS} -c $< -o $@

example1_BDADI_driver.o: example1_BDADI_driver.cpp
	${CC} ${CXXFLAGS} -c $< -o $@

example1_BDADI: example1_BDADI_driver.o argument_parser.o example1_BDADI.o reporter.o
	${CC} ${CXXFLAGS} $^ -o $@

example2.o: example2.cpp FAI_2DFSDEsolver.h
	${CC} ${CXXFLAGS} -c $< -o $@

example2_driver.o: example2_driver.cpp
	${CC} ${CXXFLAGS} -c $< -o $@

example2: example2_driver.o argument_parser.o example2.o kiss_fft.o reporter.o
	${CC} ${CXXFLAGS} $^ -o $@

example2_BFSMGM.o: example2_BFSMGM.cpp FAI_2DFSDEsolver.h
	${CC} ${CXXFLAGS} -c $< -o $@

example2_BFSMGM_driver.o: example2_BFSMGM_driver.cpp
	${CC} ${CXXFLAGS} -c $< -o $@

example2_BFSMGM: example2_BFSMGM_driver.o argument_parser.o example2_BFSMGM.o kiss_fft.o reporter.o
	${CC} ${CXXFLAGS} $^ -o $@

example3.o: example3.cpp FAI_2DFSDEsolver.h
	${CC} ${CXXFLAGS} -c $< -o $@

example3_driver.o: example3_driver.cpp
	${CC} ${CXXFLAGS} -c $< -o $@

example3: example3_driver.o argument_parser.o example3.o kiss_fft.o reporter.o
	${CC} ${CXXFLAGS} $^ -o $@

clean:
	@rm -f example1 example1_BDADI example2 example2_BFSMGM example3 *.o
