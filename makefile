CXX = g++
CXXFLAGS = -O2

create_dmats: create_dmats.cpp
	g++ -o create_dmats create_dmats.cpp

naive: naive.cpp
	$(CXX) $(CXXFLAGS) -o naive naive.cpp -g

run:
	./create_dmats
