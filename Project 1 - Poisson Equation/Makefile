INC=-Iinclude
LIB=-ljsoncpp -llapacke
OPTIMIZE=-O2 -O3 -Ofast

all: solve test

run: solve test
	bash ./test.sh

test: $(JSONLIB) matrix.o mathExpr.o test.o BVP.o test.o 
	g++ -o test matrix.o mathExpr.o solver.o BVP.o test.o $(LIB) $(OPTIMIZE)

solve: $(JSONLIB) matrix.o mathExpr.o solver.o BVP.o solve.o 
	g++ -o solve matrix.o mathExpr.o solver.o BVP.o solve.o $(LIB) $(OPTIMIZE)

matrix.o: src/matrix.cpp include/matrix.h
	g++ -c src/matrix.cpp $(INC) $(OPTIMIZE)

mathExpr.o: src/mathExpr.cpp include/mathExpr.h
	g++ -c src/mathExpr.cpp $(INC) $(OPTIMIZE)

solver.o: src/solver.cpp include/solver.h src/matrix.cpp include/matrix.h src/mathExpr.cpp include/mathExpr.h
	g++ -c src/solver.cpp $(INC) $(OPTIMIZE)

BVP.o: src/BVP.cpp include/BVP.h src/solver.cpp include/solver.h src/matrix.cpp include/matrix.h src/mathExpr.cpp include/mathExpr.h
	g++ -c src/BVP.cpp $(INC) $(OPTIMIZE)

solve.o: solve.cpp src/BVP.cpp include/BVP.h src/solver.cpp include/solver.h src/matrix.cpp include/matrix.h src/mathExpr.cpp include/mathExpr.h
	g++ -c solve.cpp $(INC) $(OPTIMIZE)

test.o: test.cpp src/BVP.cpp include/BVP.h src/solver.cpp include/solver.h src/matrix.cpp include/matrix.h src/mathExpr.cpp include/mathExpr.h
	g++ -c test.cpp $(INC) $(OPTIMIZE)

clear:
	rm solve test *.o

clean:
	if [ -f "result.txt" ];then rm -d result.txt; fi
	rm solve test *.o
	rm -r jsoncpp-1.9.5/build