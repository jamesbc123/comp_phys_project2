all: compile run

compile:
	g++ -o main.out main.cpp solver.cpp -larmadillo

run:
	./main.out



# unit_testing.cpp (perform unit tests):
all_test: compile_test link_test run_test

compile_test:
	g++ -c -Wall unit_testing.cpp unit_testing_functions.cpp solver.cpp # -c: compile only.

link_test:
	g++ -o unit_testing.out unit_testing.o unit_testing_functions.o solver.o -larmadillo

run_test:
	./unit_testing.out
