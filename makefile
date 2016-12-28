main: main.cpp main.hpp
	g++  main.cpp -o main -pthread -std=gnu++11
clean:
	$(RM) main
