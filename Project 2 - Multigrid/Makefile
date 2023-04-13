JSONLIB=jsoncpp-1.9.5/build/lib/libjsoncpp.a
INC=-Ijsoncpp-1.9.5/include/json -Iinclude
LIB=-Ljsoncpp-1.9.5/build/lib -ljsoncpp
OPTIMIZE=-O2 -O3 -Ofast

all: solve

run: solve
	bash ./test.sh

solve: $(JSONLIB) matrix.o sparseMatrix.o mathExpr.o solver.o BVP.o solve.o 
	g++ -o solve matrix.o sparseMatrix.o mathExpr.o solver.o BVP.o solve.o $(LIB) $(OPTIMIZE)

$(JSONLIB):
	cd jsoncpp-1.9.5 && mkdir build
	cd jsoncpp-1.9.5/build && cmake ..
	cd jsoncpp-1.9.5/build && make

matrix.o: src/matrix.cpp include/matrix.h
	g++ -c src/matrix.cpp $(INC) $(OPTIMIZE)

sparseMatrix.o: src/matrix.cpp include/matrix.h src/sparseMatrix.cpp include/sparseMatrix.h
	g++ -c src/sparseMatrix.cpp $(INC) $(OPTIMIZE)

solver.o: src/solver.cpp include/solver.h src/matrix.cpp include/matrix.h src/sparseMatrix.cpp include/sparseMatrix.h
	g++ -c src/solver.cpp $(INC) $(OPTIMIZE)

mathExpr.o: src/mathExpr.cpp include/mathExpr.h
	g++ -c src/mathExpr.cpp $(INC) $(OPTIMIZE)

BVP.o: src/BVP.cpp include/BVP.h src/solver.cpp include/solver.h src/matrix.cpp include/matrix.h src/sparseMatrix.cpp include/sparseMatrix.h src/mathExpr.cpp include/mathExpr.h
	g++ -c src/BVP.cpp $(INC) $(OPTIMIZE)

solve.o: solve.cpp src/BVP.cpp include/BVP.h
	g++ -c solve.cpp $(INC) $(OPTIMIZE)

clear:
	if [ -f "result.txt" ];then rm -d result.txt; fi
	rm solve *.o

clean:
	if [ -f "result.txt" ];then rm -d result.txt; fi
	rm solve *.o
	rm -r jsoncpp-1.9.5/build
