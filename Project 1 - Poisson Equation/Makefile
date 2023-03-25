JSONLIB=jsoncpp-1.9.5/build/lib/libjsoncpp.a
INC=-Ijsoncpp-1.9.5/include/json -Iinclude
LIB=-Ljsoncpp-1.9.5/build/lib -ljsoncpp
OPTIMIZE=-O2 -O3 -Ofast

all: $(JSONLIB) matrix.o mathExpr.o solver.o BVP.o test.o 
	g++ -o test matrix.o mathExpr.o solver.o BVP.o test.o $(LIB) $(OPTIMIZE)

$(JSONLIB):
	cd jsoncpp-1.9.5 && mkdir build
	cd jsoncpp-1.9.5/build && cmake ..
	cd jsoncpp-1.9.5/build && make

matrix.o: src/matrix.cpp include/matrix.h
	g++ -c src/matrix.cpp $(INC) $(OPTIMIZE)

mathExpr.o: src/mathExpr.cpp include/mathExpr.h
	g++ -c src/mathExpr.cpp $(INC) $(OPTIMIZE)

solver.o: src/solver.cpp include/solver.h src/matrix.cpp include/matrix.h src/mathExpr.cpp include/mathExpr.h
	g++ -c src/solver.cpp $(INC) $(OPTIMIZE)

BVP.o: src/BVP.cpp include/BVP.h src/solver.cpp include/solver.h src/matrix.cpp include/matrix.h src/mathExpr.cpp include/mathExpr.h
	g++ -c src/BVP.cpp $(INC) $(OPTIMIZE)

test.o: test.cpp src/BVP.cpp include/BVP.h src/solver.cpp include/solver.h src/matrix.cpp include/matrix.h src/mathExpr.cpp include/mathExpr.h
	g++ -c test.cpp $(INC) $(OPTIMIZE)

clear:
	rm test *.o

clean:
	if [ -f "result.txt" ];then rm -d result.txt; fi
	rm test *.o
	rm -r jsoncpp-1.9.5/build