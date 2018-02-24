CFLAGS=-O3 -D_FORCE_INLINES -std=c++11 -D_MWAITXINTRIN_H_INCLUDED


#pagerank_test: pagerank_test.cpp table.cpp pagerank.cpp table.h
#	g++ $(CFLAGS) -o pagerank_test pagerank_test.cpp table.cpp

pagerank_test: pagerank_test.cpp table.cpp pagerank.cpp table.cuh
	/vol/cuda/9.0.103-RC1/bin/nvcc $(CFLAGS) -o pagerank_test pagerank_test.cpp table.cu 

pagerank: pagerank.cpp table.cpp table.cuh
	nvcc $(CFLAGS) --compiler-options -Wall -o pagerank pagerank.cpp table.cu 

all-tests: all-tests.txt pagerank_test
	./pagerank_test all-tests.txt

small-test: small pagerank_test
	./pagerank_test small

medium-test: medium pagerank_test
	./pagerank_test medium

large-test: large pagerank_test
	./pagerank_test large

enormous-test: enormous pagerank_test
	./pagerank_test enormous

ginormous-test: ginormous pagerank_test
	./pagerank_test ginormous

bull-test: bull pagerank_test
	./pagerank_test bull
clean:
	rm -f pagerank pagerank_test
