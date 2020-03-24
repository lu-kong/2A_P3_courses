CC = g++

first: first.cpp
		$(CC) -o first first.cpp
clean: 
		rm -f first *~		