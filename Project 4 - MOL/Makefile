INC=-Iinclude
LIB=-ljsoncpp
OPTIMIZE=-O2 -O3 -Ofast

story: heat advection trueSol
	cd plotter && bash picgen.sh
	cd report && make

trueSol: src/trueSol.cpp
	g++ src/trueSol.cpp -o trueSol $(OPTIMIZE)

heat: heat.o IVP.o LMMs.o RungeKutta.o matrix.o NonLinearSolver.o Polynomial.o
	g++ heat.o IVP.o LMMs.o RungeKutta.o matrix.o NonLinearSolver.o Polynomial.o $(LIB) -o heat $(OPTIMIZE)

advection: advection.o matrix.o
	g++ advection.o matrix.o $(LIB) -o advection $(OPTIMIZE)

$(JSONLIB):
	cd jsoncpp-1.9.5 && mkdir build
	cd jsoncpp-1.9.5/build && cmake ..
	cd jsoncpp-1.9.5/build && make

heat.o: src/heat.cpp
	g++ -c src/heat.cpp -o heat.o $(INC) $(OPTIMIZE)

advection.o: src/advection.cpp
	g++ -c src/advection.cpp -o advection.o $(INC) $(OPTIMIZE)

solve.o: solve.cpp functions/*.h
	g++ -c solve.cpp $(INC) -Ifunctions $(OPTIMIZE)

matrix.o: src/matrix.cpp include/matrix.h
	g++ -c src/matrix.cpp $(INC) $(OPTIMIZE)

IVP.o: src/IVP.cpp include/IVP.h src/matrix.cpp include/matrix.h src/Polynomial.cpp include/Polynomial.h
	g++ -c src/IVP.cpp $(INC) $(OPTIMIZE)

NonLinearSolver.o: src/NonLinearSolver.cpp include/NonLinearSolver.h src/matrix.cpp include/matrix.h
	g++ -c src/NonLinearSolver.cpp $(INC) $(OPTIMIZE)

RungeKutta.o: src/RungeKutta.cpp include/RungeKutta.h src/IVP.cpp include/IVP.h src/NonLinearSolver.cpp include/NonLinearSolver.h src/Polynomial.cpp include/Polynomial.h
	g++ -c src/RungeKutta.cpp $(INC) $(OPTIMIZE)

Polynomial.o: src/Polynomial.cpp include/Polynomial.h src/matrix.cpp include/matrix.h
	g++ -c src/Polynomial.cpp $(INC) $(OPTIMIZE)

LMMs.o: src/LMMs.cpp include/LMMs.h src/IVP.cpp include/IVP.h src/NonLinearSolver.cpp include/NonLinearSolver.h
	g++ -c src/LMMs.cpp $(INC) $(OPTIMIZE)

clean:
	rm *.o heat advection trueSol
	cd report && make clean
