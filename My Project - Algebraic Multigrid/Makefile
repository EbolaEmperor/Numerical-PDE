all: test

test: test.o amgSolver.o sparseMatrix.o matrix.o avl.o
	g++ test.o sparseMatrix.o matrix.o avl.o amgSolver.o -o test -O2 -O3 -Ofast

test.o: test.cpp
	g++ -c test.cpp -Iinclude -O2 -O3 -Ofast

avl.o: src/avl.cpp include/avl.h
	g++ -c src/avl.cpp -Iinclude -O2 -O3 -Ofast

matrix.o: src/matrix.cpp include/matrix.h
	g++ -c src/matrix.cpp -Iinclude -O2 -O3 -Ofast

sparseMatrix.o: src/sparseMatrix.cpp include/sparseMatrix.h src/matrix.cpp include/matrix.h
	g++ -c src/sparseMatrix.cpp -Iinclude -O2 -O3 -Ofast

amgSolver.o: src/amgSolver.cpp include/amgSolver.h src/avl.cpp include/avl.h src/sparseMatrix.cpp include/sparseMatrix.h src/matrix.cpp include/matrix.h
	g++ -c src/amgSolver.cpp -Iinclude -O2 -O3 -Ofast

clean:
	rm *.o test