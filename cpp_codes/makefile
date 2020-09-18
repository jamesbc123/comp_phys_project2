all: compile run

compile:
	g++ -o main.out main.cpp solver.cpp

run:
	./main.out


all_test: compile_test link_test run_test

#compile_test:
#	g++ -o unit_testing.out unit_testing.cpp

#run_test:
#	./unit_testing.out


compile_test:
	g++ -c -Wall unit_testing.cpp # -c: compile only.

link_test:
	g++ -o unit_testing.out unit_testing.o -larmadillo

run_test:
	./unit_testing.out
