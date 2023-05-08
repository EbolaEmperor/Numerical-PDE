all: test

run: test
	bash ./run.sh

test: test.o IVP.o LMMs.o RungeKutta.o matrix.o NonLinearSolver.o Polynomial.o
	g++ test.o IVP.o LMMs.o RungeKutta.o matrix.o NonLinearSolver.o Polynomial.o -o test -O2 -O3 -Ofast

test.o: test.cpp
	g++ -c test.cpp -Iinclude -O2 -O3 -Ofast

matrix.o: src/matrix.cpp include/matrix.h
	g++ -c src/matrix.cpp -Iinclude -O2 -O3 -Ofast

IVP.o: src/IVP.cpp include/IVP.h src/matrix.cpp include/matrix.h src/Polynomial.cpp include/Polynomial.h
	g++ -c src/IVP.cpp -Iinclude -O2 -O3 -Ofast

NonLinearSolver.o: src/NonLinearSolver.cpp include/NonLinearSolver.h src/matrix.cpp include/matrix.h
	g++ -c src/NonLinearSolver.cpp -Iinclude -O2 -O3 -Ofast

RungeKutta.o: src/RungeKutta.cpp include/RungeKutta.h src/IVP.cpp include/IVP.h src/NonLinearSolver.cpp include/NonLinearSolver.h src/Polynomial.cpp include/Polynomial.h
	g++ -c src/RungeKutta.cpp -Iinclude -O2 -O3 -Ofast

Polynomial.o: src/Polynomial.cpp include/Polynomial.h src/matrix.cpp include/matrix.h
	g++ -c src/Polynomial.cpp -Iinclude -O2 -O3 -Ofast

LMMs.o: src/LMMs.cpp include/LMMs.h src/IVP.cpp include/IVP.h src/NonLinearSolver.cpp include/NonLinearSolver.h
	g++ -c src/LMMs.cpp -Iinclude -O2 -O3 -Ofast

clean:
	$(shell if [ -e Backward\ Eular ];then rm Backward\ Eular; fi)
	$(shell if [ -e Trapezoidal ];then rm Trapezoidal; fi)
	rm test *.o
