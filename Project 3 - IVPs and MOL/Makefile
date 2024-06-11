INC=-Iinclude
LIB=-ljsoncpp
OPTIMIZE=-O2 -O3 -Ofast

story: solve
	bash picgen.sh
	cd report && make
	cd report && make release

run: solve
	./solve tests/*.json

solve: solve.o IVP.o LMMs.o RungeKutta.o SymplecticMethods.o matrix.o NonLinearSolver.o Polynomial.o functionFactory.o
	g++ solve.o IVP.o LMMs.o RungeKutta.o SymplecticMethods.o matrix.o NonLinearSolver.o Polynomial.o functionFactory.o $(LIB) -o solve $(OPTIMIZE)

heat: MOL/heat.o IVP.o LMMs.o RungeKutta.o matrix.o NonLinearSolver.o Polynomial.o
	g++ MOL/heat.o IVP.o LMMs.o RungeKutta.o matrix.o NonLinearSolver.o Polynomial.o $(LIB) -o MOL/heat $(OPTIMIZE)

advection: MOL/advection.o matrix.o
	g++ MOL/advection.o matrix.o $(LIB) -o MOL/advection $(OPTIMIZE)

MOL/heat.o: MOL/heat.cpp
	g++ -c MOL/heat.cpp -o MOL/heat.o $(INC) $(OPTIMIZE)

MOL/advection.o: MOL/advection.cpp
	g++ -c MOL/advection.cpp -o MOL/advection.o $(INC) $(OPTIMIZE)

solve.o: solve.cpp functions/*.h
	g++ -c solve.cpp $(INC) -Ifunctions $(OPTIMIZE)

functionFactory.o: src/functionFactory.cpp include/functionFactory.h
	g++ -c src/functionFactory.cpp $(INC) $(OPTIMIZE)

matrix.o: src/matrix.cpp include/matrix.h
	g++ -c src/matrix.cpp $(INC) $(OPTIMIZE)

IVP.o: src/IVP.cpp include/IVP.h src/matrix.cpp include/matrix.h src/Polynomial.cpp include/Polynomial.h
	g++ -c src/IVP.cpp $(INC) $(OPTIMIZE)

NonLinearSolver.o: src/NonLinearSolver.cpp include/NonLinearSolver.h src/matrix.cpp include/matrix.h
	g++ -c src/NonLinearSolver.cpp $(INC) $(OPTIMIZE)

RungeKutta.o: src/RungeKutta.cpp include/RungeKutta.h src/IVP.cpp include/IVP.h src/NonLinearSolver.cpp include/NonLinearSolver.h src/Polynomial.cpp include/Polynomial.h
	g++ -c src/RungeKutta.cpp $(INC) $(OPTIMIZE)

SymplecticMethods.o: src/SymplecticMethods.cpp include/SymplecticMethods.h src/IVP.cpp include/IVP.h
	g++ -c src/SymplecticMethods.cpp $(INC) $(OPTIMIZE)

Polynomial.o: src/Polynomial.cpp include/Polynomial.h src/matrix.cpp include/matrix.h
	g++ -c src/Polynomial.cpp $(INC) $(OPTIMIZE)

LMMs.o: src/LMMs.cpp include/LMMs.h src/IVP.cpp include/IVP.h src/NonLinearSolver.cpp include/NonLinearSolver.h
	g++ -c src/LMMs.cpp $(INC) $(OPTIMIZE)

clean:
	$(shell if [ -e solve ];then rm solve; fi)
	$(shell if [ -e MOL/heat ];then rm MOL/heat MOL/*.txt MOL/*.o; fi)
	$(shell if [ -e MOL/advection ];then rm MOL/advection; fi)
	$(shell if [ -e result.txt ];then rm *.txt; fi)
	rm *.o
	cd report && make clean
	cd report && cd figures && mv 11-1.eps 11-1.epsback && rm *.eps 5-5.png && mv 11-1.epsback 11-1.eps
