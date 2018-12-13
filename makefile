CXX = g++
CXXFLAGS = -O2

create_dmats: create_dmats.cpp
	g++ -o create_dmats create_dmats.cpp

naive: naive.cpp
	$(CXX) $(CXXFLAGS) -o naive naive.cpp

slow: naive_v2.cpp
	$(CXX) -o slow naive_v2.cpp

opt1: naive_v2.cpp
	$(CXX) -O1 -o opt1 naive_v2.cpp

opt2: naive_v2.cpp
	$(CXX) -O2 -o opt2 naive_v2.cpp

opt3: naive_v2.cpp
	$(CXX) -O3 -o opt3 naive_v2.cpp

wall: parallel_cpu_wall.cu
	nvcc parallel_cpu_wall.cu -o wall -arch=sm_35 -D_FORCE_INLINES

parallel: parallel.cu
	nvcc parallel.cu -o parallel -arch=sm_35 -D_FORCE_INLINES

run:
	./create_dmats
