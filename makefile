CXX = g++
CXXFLAGS = -O2

create_dmats: create_dmats.cpp
	g++ -o create_dmats create_dmats.cpp

naive: naive.cpp
	$(CXX) $(CXXFLAGS) -o naive naive.cpp

slow: naive.cpp
	$(CXX) -o slow naive.cpp

opt1: naive.cpp
	$(CXX) -O1 -o opt1 naive.cpp

opt2: naive.cpp
	$(CXX) -O2 -o opt2 naive.cpp

opt3: naive.cpp
	$(CXX) -O3 -o opt3 naive.cpp

v2: naive_v2.cpp
	$(CXX) -O3 -o naive2 naive_v2.cpp

parallel: parallel.cu
	nvcc parallel.cu -o parallel -arch=sm_35 -D_FORCE_INLINES

run:
	./create_dmats
