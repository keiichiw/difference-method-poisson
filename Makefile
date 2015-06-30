CC=g++
FLG=-Wall -Wextra -std=c++11
LIB=~/lib/eigen/

all:
	$(CC) $(FLG) -I $(LIB) main.cpp -o main

clean:
	rm -f main *.o *~
