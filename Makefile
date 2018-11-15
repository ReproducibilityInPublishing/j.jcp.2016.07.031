all: example1 example1_BDADI example2 example2_BFSMGM example3

ifndef CXX
CXX := g++
endif

ifndef CC
CC := gcc
endif

CXXFLAGS := -g -O3 --std=c++11

argument_parser.o: argument_parser.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

kiss_fft.o: kiss_fft.c
	${CC} ${CXXFLAGS} -c $< -o $@

reporter.o: reporter.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

example1.o: example1.cpp FAI_2DFSDEsolver.h
	${CXX} ${CXXFLAGS} -c $< -o $@

example1_driver.o: example1_driver.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

example1: example1_driver.o argument_parser.o example1.o kiss_fft.o reporter.o
	${CXX} ${CXXFLAGS} $^ -o $@

example1_BDADI.o: example1_BDADI.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

example1_BDADI_driver.o: example1_BDADI_driver.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

example1_BDADI: example1_BDADI_driver.o argument_parser.o example1_BDADI.o reporter.o
	${CXX} ${CXXFLAGS} $^ -o $@

example2.o: example2.cpp FAI_2DFSDEsolver.h
	${CXX} ${CXXFLAGS} -c $< -o $@

example2_driver.o: example2_driver.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

example2: example2_driver.o argument_parser.o example2.o kiss_fft.o reporter.o
	${CXX} ${CXXFLAGS} $^ -o $@

example2_BFSMGM.o: example2_BFSMGM.cpp FAI_2DFSDEsolver.h
	${CXX} ${CXXFLAGS} -c $< -o $@

example2_BFSMGM_driver.o: example2_BFSMGM_driver.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

example2_BFSMGM: example2_BFSMGM_driver.o argument_parser.o example2_BFSMGM.o kiss_fft.o reporter.o
	${CXX} ${CXXFLAGS} $^ -o $@

example3.o: example3.cpp FAI_2DFSDEsolver.h
	${CXX} ${CXXFLAGS} -c $< -o $@

example3_driver.o: example3_driver.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

example3: example3_driver.o argument_parser.o example3.o kiss_fft.o reporter.o
	${CXX} ${CXXFLAGS} $^ -o $@

clean:
	@rm -f example1 example1_BDADI example2 example2_BFSMGM example3 *.o
