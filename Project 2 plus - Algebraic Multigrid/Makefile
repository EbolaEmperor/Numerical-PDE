all: test1 test2 test3

test3: test3.o amgSolver.o sparseMatrix.o matrix.o avl.o
	g++ test3.o sparseMatrix.o matrix.o avl.o amgSolver.o -o test3 -O2 -O3 -Ofast

test3.o: test3.cpp
	g++ -c test3.cpp -Iinclude -O2 -O3 -Ofast

test2: test2.o amgSolver.o sparseMatrix.o matrix.o avl.o
	g++ test2.o sparseMatrix.o matrix.o avl.o amgSolver.o -o test2 -O2 -O3 -Ofast

test2.o: test2.cpp
	g++ -c test2.cpp -Iinclude -O2 -O3 -Ofast

test1: test1.o amgSolver.o sparseMatrix.o matrix.o avl.o
	g++ test1.o sparseMatrix.o matrix.o avl.o amgSolver.o -o test1 -O2 -O3 -Ofast

test1.o: test1.cpp
	g++ -c test1.cpp -Iinclude -O2 -O3 -Ofast

avl.o: src/avl.cpp include/avl.h
	g++ -c src/avl.cpp -Iinclude -O2 -O3 -Ofast

matrix.o: src/matrix.cpp include/matrix.h
	g++ -c src/matrix.cpp -Iinclude -O2 -O3 -Ofast

sparseMatrix.o: src/sparseMatrix.cpp include/sparseMatrix.h src/matrix.cpp include/matrix.h
	g++ -c src/sparseMatrix.cpp -Iinclude -O2 -O3 -Ofast

amgSolver.o: src/amgSolver.cpp include/amgSolver.h src/avl.cpp include/avl.h src/sparseMatrix.cpp include/sparseMatrix.h src/matrix.cpp include/matrix.h
	g++ -c src/amgSolver.cpp -Iinclude -O2 -O3 -Ofast

clean:
	rm *.o test1 test2 test3
