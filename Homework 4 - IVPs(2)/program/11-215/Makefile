all: test

run: test
	./test

test: test.cpp polynomial.h fraction.h bigint.h
	g++ test.cpp -o test -O2 -O3 -Ofast

clean:
	rm test